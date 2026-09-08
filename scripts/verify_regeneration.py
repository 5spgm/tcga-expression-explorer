"""2つの出力ディレクトリの遺伝子JSONを全遺伝子で突き合わせる。

  python3 verify_regeneration.py <公開側ディレクトリ> <再生成側ディレクトリ> [出力JSON]

同じ入力・同じ処理なら完全一致するはず(構造の取り違えやビット化けを検出する)。
入力を変えた世代については相対差・絶対差の最大値を測る。
BRCA中期のように Excel往復の丸めが消える場合、差は倍精度の丸め幅に収まる。
"""
import json, os, sys, math
from multiprocessing import Pool

A_DIR, B_DIR = sys.argv[1], sys.argv[2]
OUT_JSON = sys.argv[3] if len(sys.argv) > 3 else "verify_regeneration.json"


def walk(obj, path=""):
    """入れ子の dict/list から (パス, 数値リスト) を取り出す"""
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from walk(v, f"{path}/{k}")
    elif isinstance(obj, list):
        yield path, obj


def compare(fname):
    pa, pb = os.path.join(A_DIR, fname), os.path.join(B_DIR, fname)
    if not os.path.exists(pb):
        return {"missing": fname}
    try:
        with open(pa, encoding="utf-8") as fh:
            a = json.load(fh)
    except Exception as e:  # noqa: BLE001
        return {"unreadable": f"A/{fname}: {type(e).__name__} {e}"}
    try:
        with open(pb, encoding="utf-8") as fh:
            b = json.load(fh)
    except Exception as e:  # noqa: BLE001
        return {"unreadable": f"B/{fname}: {type(e).__name__} {e}"}

    res = {"file": fname, "meta_ok": all(a.get(k) == b.get(k) for k in
                                         ("gene_id", "gene_symbol", "gene_type")),
           "gen_keys_ok": set(a) == set(b),
           "exact_gen": [], "diff_gen": [], "len_mismatch": [],
           "max_rel": 0.0, "max_abs": 0.0, "n_values": 0}

    for gen in ("new", "mid", "old"):
        if gen not in a and gen not in b:
            # その世代を持たない遺伝子(旧世代は約4万遺伝子で欠ける)。両側で欠けて
            # いるのは一致であって不一致ではない。
            continue
        if gen not in a or gen not in b:
            res["gen_keys_ok"] = False
            continue
        la, lb = dict(walk(a[gen])), dict(walk(b[gen]))
        if set(la) != set(lb):
            res["gen_keys_ok"] = False
            continue
        identical = True
        for k in la:
            va, vb = la[k], lb[k]
            if len(va) != len(vb):
                res["len_mismatch"].append(f"{gen}{k}: {len(va)} vs {len(vb)}")
                identical = False
                continue
            for x, y in zip(va, vb):
                if x is None or y is None:
                    if x is not y:
                        identical = False
                    continue
                res["n_values"] += 1
                if x != y:
                    identical = False
                    d = abs(x - y)
                    res["max_abs"] = max(res["max_abs"], d)
                    denom = max(abs(x), abs(y))
                    if denom > 0:
                        res["max_rel"] = max(res["max_rel"], d / denom)
                    # new / old は同じ入力・同じ丸めなので一致するはず。
                    # 不一致は値そのものを残して後で原因を判断する。
                    if gen in ("new", "old"):
                        res.setdefault("hard_diff", []).append(
                            f"{fname} {gen}{k} 現行={x!r} 再生成={y!r}")
        (res["exact_gen"] if identical else res["diff_gen"]).append(gen)
    return res


def main():
    files = sorted(f for f in os.listdir(A_DIR)
                   if f.endswith(".json") and not f.startswith("_"))
    print(f"対象 {len(files):,} 遺伝子JSON", flush=True)
    agg = {"n": 0, "missing": [], "meta_bad": [], "keys_bad": [], "len_bad": [],
           "unreadable": [], "hard_diff": [],
           "max_rel": 0.0, "max_abs": 0.0, "n_values": 0,
           "exact": {"new": 0, "mid": 0, "old": 0},
           "diff": {"new": 0, "mid": 0, "old": 0},
           "worst_rel": None}
    with Pool(16) as pool:
        for r in pool.imap_unordered(compare, files, chunksize=64):
            if "missing" in r:
                agg["missing"].append(r["missing"]); continue
            if "unreadable" in r:
                agg["unreadable"].append(r["unreadable"]); continue
            agg["hard_diff"] += r.get("hard_diff", [])
            agg["n"] += 1
            agg["n_values"] += r["n_values"]
            if not r["meta_ok"]:
                agg["meta_bad"].append(r["file"])
            if not r["gen_keys_ok"]:
                agg["keys_bad"].append(r["file"])
            agg["len_bad"] += [f"{r['file']} {m}" for m in r["len_mismatch"]]
            for g in r["exact_gen"]:
                agg["exact"][g] += 1
            for g in r["diff_gen"]:
                agg["diff"][g] += 1
            if r["max_rel"] > agg["max_rel"]:
                agg["max_rel"] = r["max_rel"]
                agg["worst_rel"] = r["file"]
            agg["max_abs"] = max(agg["max_abs"], r["max_abs"])
            if agg["n"] % 10000 == 0:
                print(f"  {agg['n']:,} 件", flush=True)

    print(f"\n比較した遺伝子: {agg['n']:,}  値の総数: {agg['n_values']:,}")
    print(f"再生成側に無い: {len(agg['missing'])}  メタ情報不一致: {len(agg['meta_bad'])}  "
          f"構造不一致: {len(agg['keys_bad'])}  配列長不一致: {len(agg['len_bad'])}")
    for g in ("new", "mid", "old"):
        print(f"  {g:4} 完全一致 {agg['exact'][g]:,} 遺伝子 / 差あり {agg['diff'][g]:,} 遺伝子")
    print(f"最大相対差: {agg['max_rel']:.3e}  (遺伝子 {agg['worst_rel']})")
    print(f"最大絶対差: {agg['max_abs']:.3e}")
    print(f"読めなかったファイル: {len(agg['unreadable'])}")
    for u in agg["unreadable"][:5]:
        print("   ", u)
    print(f"new/old の値の不一致: {len(agg['hard_diff'])} 箇所")
    for h in agg["hard_diff"][:10]:
        print("   ", h)
    for label in ("missing", "meta_bad", "keys_bad", "len_bad"):
        if agg[label]:
            print(f"  [{label}] 例: {agg[label][:5]}")
    json.dump({k: v for k, v in agg.items() if k not in ("missing", "meta_bad", "keys_bad", "len_bad", "unreadable", "hard_diff")}
              | {k: agg[k][:50] for k in ("missing", "meta_bad", "keys_bad", "len_bad",
                                          "unreadable", "hard_diff")},
              open(OUT_JSON, "w"), indent=1)
    print(f"→ {OUT_JSON}")


if __name__ == "__main__":
    main()
