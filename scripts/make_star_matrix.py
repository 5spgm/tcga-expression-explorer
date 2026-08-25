#!/usr/bin/env python3
"""
make_star_matrix.py
===================
GDCからダウンロードした新世代の STAR-Counts ファイル
(`*.rna_seq.augmented_star_gene_counts.tsv`)を1つの発現行列にまとめる。
`fpkmlist.R` / `20250913_RNAseq_list.R` の置き換え。

出力は `gene_id, gene_name, gene_type, TCGA-XX-XXXX-01A, ...` の形のTSVで、
preprocess_tcga.py の `--new-tpm` などにそのまま渡せる。

## 使い方

    python3 make_star_matrix.py \\
        --sample-sheet gdc_sample_sheet.2023-02-21.tsv \\
        --data-dir     ./gdc_download \\
        --value        tpm_unstranded \\
        --output       TCGA_PAAD_TPM.tsv

1ファイルを読むだけで9列すべてが手に入るので、複数の値を一度に出せる:

    --value tpm_unstranded,fpkm_unstranded,fpkm_uq_unstranded
    -> TCGA_PAAD_TPM.tsv / _FPKM.tsv / _FPKM_UQ.tsv を同時に生成

## 元のRスクリプトから変えた点

1. **merge() の繰り返しをやめた**
   `for(i in 2:nrow) { fpkm1 <- merge(fpkm1, fpkm) }` は毎回テーブル全体を
   作り直すためO(n^2)になり、1,000検体規模で数時間かかる。
   STARの出力は全ファイルで遺伝子の並びが同一なので、最初の1ファイルの
   並びを基準にして、以降は値の列だけを事前確保した配列へ流し込む。
   並びが違うファイルがあれば検出して停止する(黙って壊れるより安全)。

2. **検体を落とさない**
   `subset(manifest, [,8]=="Primary Tumor")` は Metastatic / Recurrent を
   警告なく捨てていた。ここでは全てのsample typeを残し、原発巣か転移巣かの
   判定は preprocess_tcga.py 側に任せる(そちらで既定除外＋チェックボックス
   切替ができる)。

3. **`grep("-01A", ...)` を使わない**
   あれは「複数ファイルを持つ患者で01Aが無い場合、その患者ごと消える」
   という副作用があった。代わりに、同一患者×同一sample typeに複数vialが
   ある場合だけ1つに絞る(既定はvial文字が若い方)。

4. **患者単位ではなく検体単位で重複除去**
   `!duplicated(man3[,6])` は6列目=Case IDでの除去になっていた。

5. **biomaRt を使わない**
   STARの出力には `gene_name` 列が入っているので、Ensembl ID から
   シンボルを引き直す必要がない。ネット接続もbiomaRtの障害も無関係になる。

6. **`.xls` ではなくTSVで書く**
   元のスクリプトは write.table でTSVを書いて拡張子だけ .xls にしていた。
   ここでは素直に .tsv とし、後段の読み込みも速くなる。
"""

import argparse
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

VALUE_COLUMNS = [
    "unstranded", "stranded_first", "stranded_second",
    "tpm_unstranded", "fpkm_unstranded", "fpkm_uq_unstranded",
]
# --value に対応する出力ファイル名の接尾辞
SUFFIX = {
    "tpm_unstranded": "TPM",
    "fpkm_unstranded": "FPKM",
    "fpkm_uq_unstranded": "FPKM_UQ",
    "unstranded": "counts",
    "stranded_first": "counts_first",
    "stranded_second": "counts_second",
}

BARCODE_TAIL = "-"  # Sample ID は TCGA-XX-XXXX-01A の形


def sample_type_code(sample_id: str) -> str:
    """TCGA-XX-XXXX-01A -> '01'"""
    parts = str(sample_id).split("-")
    return parts[3][:2] if len(parts) >= 4 else ""


def vial_letter(sample_id: str) -> str:
    parts = str(sample_id).split("-")
    return parts[3][2:3] if len(parts) >= 4 and len(parts[3]) > 2 else ""


def patient_id(sample_id: str) -> str:
    return str(sample_id)[:12]


def select_files(sheet: pd.DataFrame, prefer_vial: str) -> pd.DataFrame:
    """1検体につき1ファイルを選び、同一患者×同一sample typeの重複vialを
    1つに絞る。sample typeによる取捨選択はここでは行わない。"""
    print(f"  sample sheet: {len(sheet)} 行 / "
          f"{sheet['Case ID'].nunique()} 患者 / {sheet['Sample ID'].nunique()} 検体")
    print("  sample type の内訳:")
    for val, n in sheet["Sample Type"].value_counts().items():
        print(f"      {val}: {n}")

    # (1) 同一 Sample ID に複数ファイルがある場合は先頭を採用
    dup_samples = sheet["Sample ID"].value_counts()
    n_dup = int((dup_samples > 1).sum())
    picked = sheet.drop_duplicates(subset=["Sample ID"], keep="first").copy()
    if n_dup:
        print(f"  同一検体に複数ファイルがある {n_dup} 検体は先頭の1ファイルを採用しました")

    # (2) 同一患者×同一sample typeで複数vial -> 1つに絞る
    groups = defaultdict(list)
    for idx, row in picked.iterrows():
        groups[(patient_id(row["Sample ID"]), sample_type_code(row["Sample ID"]))].append(idx)

    drop_idx, notes = [], []
    for (pid, code), idxs in groups.items():
        if len(idxs) == 1:
            continue
        ordered = sorted(idxs, key=lambda i: vial_letter(picked.loc[i, "Sample ID"]))
        keep = ordered[0] if prefer_vial == "first" else ordered[-1]
        for i in ordered:
            if i != keep:
                drop_idx.append(i)
        notes.append(f"{pid}-{code}: "
                     f"{'/'.join(picked.loc[i, 'Sample ID'][-3:] for i in ordered)} -> "
                     f"{picked.loc[keep, 'Sample ID'][-3:]} を採用")
    if drop_idx:
        print(f"  同一患者×同一sample typeの重複vialを {len(notes)} 組検出しました:")
        for line in notes[:10]:
            print(f"      {line}")
        if len(notes) > 10:
            print(f"      ... 他 {len(notes) - 10} 組")
        picked = picked.drop(index=drop_idx)

    print(f"  -> 採用: {len(picked)} 検体")
    return picked


# GDCの標準ファイル名は
#   <UUID>.rna_seq.augmented_star_gene_counts.tsv
# だが、手元で短くリネームしてあるものが混在することがある。
# sample sheet 側は元の名前のままなので、UUID部分で突き合わせられるようにする。
KNOWN_SUFFIXES = [
    ".rna_seq.augmented_star_gene_counts",
    ".augmented_star_gene_counts",
    ".rna_seq",
]
SCAN_PATTERNS = ["*.tsv", "*.txt", "*.tsv.gz", "*.txt.gz"]


def file_uuid(filename: str) -> str:
    """ファイル名の先頭(最初のドットまで)をUUIDとみなす。
    リネームの有無に関わらず変わらない部分。"""
    return str(filename).split(".")[0].strip().lower()


def index_files(data_dir: Path) -> tuple[dict, dict]:
    """data-dir 以下を再帰的に走査し、
      (ファイル名 -> 実パス, UUID -> 実パス)
    の2つの対応表を作る。GDCクライアントは <File ID>/<File Name> の形で
    展開するため、サブディレクトリを掘って探す必要がある。"""
    by_name, by_uuid = {}, {}
    collisions = []
    n_files = 0
    for pattern in SCAN_PATTERNS:
        for p in data_dir.rglob(pattern):
            if not p.is_file():
                continue
            n_files += 1
            by_name.setdefault(p.name, p)
            u = file_uuid(p.name)
            if u in by_uuid and by_uuid[u] != p:
                collisions.append((u, by_uuid[u], p))
            else:
                by_uuid.setdefault(u, p)
    print(f"  {data_dir} 以下で {n_files} 件のファイルが見つかりました "
          f"(UUIDで一意に識別できたもの: {len(by_uuid)} 件)")
    if collisions:
        print(f"  [!] 同じUUIDのファイルが複数あります({len(collisions)} 件)。"
              f"先に見つかった方を使います:")
        for u, a, b in collisions[:5]:
            print(f"      {u}: {a.name} / {b.name}")
    return by_name, by_uuid


def resolve_path(wanted: str, by_name: dict, by_uuid: dict):
    """sample sheet の File Name から実ファイルを探す。
    1) 名前がそのまま一致
    2) 既知の接尾辞を外した名前で一致
    3) UUID(先頭部分)で一致  ← リネーム済みのファイルはここで拾う
    戻り値: (パス or None, 一致方法)
    """
    if wanted in by_name:
        return by_name[wanted], "name"

    stem, ext = wanted, ""
    for e in (".tsv.gz", ".txt.gz", ".tsv", ".txt"):
        if stem.endswith(e):
            stem, ext = stem[: -len(e)], e
            break
    for suf in KNOWN_SUFFIXES:
        if stem.endswith(suf):
            for e2 in (ext, ".tsv", ".txt", ".tsv.gz", ".txt.gz"):
                cand = stem[: -len(suf)] + e2
                if cand in by_name:
                    return by_name[cand], "suffix"
            break

    u = file_uuid(wanted)
    if u in by_uuid:
        return by_uuid[u], "uuid"
    return None, "missing"


def read_star_file(path: Path) -> pd.DataFrame:
    """STARの出力を読む。1行目は '# gene-model: ...' のコメント、
    先頭数行の N_unmapped 等の集計行は落とす。"""
    df = pd.read_csv(path, sep="\t", skiprows=1, dtype={"gene_id": str,
                                                        "gene_name": str,
                                                        "gene_type": str})
    df = df[~df["gene_id"].astype(str).str.startswith("N_")].reset_index(drop=True)
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-sheet", required=True, help="GDC sample sheet (TSV)")
    ap.add_argument("--data-dir", required=True,
                    help="ダウンロードしたファイルがあるディレクトリ(再帰的に探索)")
    ap.add_argument("--value", default="tpm_unstranded",
                    help=f"取り出す列。カンマ区切りで複数可。候補: {', '.join(VALUE_COLUMNS)}")
    ap.add_argument("--output", help="出力先。--value が複数のときは接尾辞が付く")
    ap.add_argument("--out-prefix", help="--output の代わりに接頭辞で指定 (例: TCGA_PAAD)")
    ap.add_argument("--prefer-vial", default="first", choices=["first", "last"],
                    help="重複vialのどちらを採るか(既定: first = A)")
    ap.add_argument("--manifest", help="gdc_manifest.txt を渡すと、"
                                       "sample sheetとの件数の食い違いを報告する")
    args = ap.parse_args()

    values = [v.strip() for v in args.value.split(",") if v.strip()]
    for v in values:
        if v not in VALUE_COLUMNS:
            sys.exit(f"--value '{v}' は不正です。候補: {', '.join(VALUE_COLUMNS)}")

    print("=" * 60)
    print("1. sample sheet から読むファイルを決める")
    print("=" * 60)
    sheet = pd.read_csv(args.sample_sheet, sep="\t", dtype=str).fillna("")
    for col in ("File Name", "Case ID", "Sample ID", "Sample Type"):
        if col not in sheet.columns:
            sys.exit(f"sample sheet に '{col}' 列がありません。列名: {list(sheet.columns)}")
    picked = select_files(sheet, args.prefer_vial)

    if args.manifest:
        man = pd.read_csv(args.manifest, sep="\t", dtype=str)
        print(f"  manifest: {len(man)} 行 / sample sheet: {len(sheet)} 行"
              f"{'  (一致)' if len(man) == len(sheet) else '  ← 食い違いあり。要確認'}")

    print()
    print("=" * 60)
    print("2. ファイルを読み込む")
    print("=" * 60)
    data_dir = Path(args.data_dir)
    if not data_dir.exists():
        sys.exit(f"--data-dir のディレクトリが存在しません: {data_dir.resolve()}\n"
                 f"  今いる場所: {Path.cwd()}\n"
                 f"  ヒント: find ~ -name '*star_gene_counts*' -maxdepth 4 | head\n"
                 f"          で実際の場所を探し、そのファイルを含むディレクトリを指定してください。")
    if not data_dir.is_dir():
        sys.exit(f"--data-dir がディレクトリではありません: {data_dir.resolve()}")

    by_name, by_uuid = index_files(data_dir)
    if not by_name:
        # 何が入っているのかを見せて、指定先の取り違えに気づけるようにする
        entries = sorted(p.name for p in data_dir.iterdir())[:10]
        n_tar = len(list(data_dir.rglob("*.tar.gz")))
        print(f"\n  {data_dir.resolve()} の中身(先頭10件):")
        for e in entries:
            print(f"      {e}")
        if not entries:
            print("      (空です)")
        if n_tar:
            print(f"\n  [!] .tar.gz が {n_tar} 件あります。展開が必要かもしれません:")
            print(f"      cd {data_dir} && tar xzf *.tar.gz")
        sys.exit("\n読み込めるファイル(.tsv / .txt)が1件もありませんでした。"
                 "--data-dir の指定先を確認してください。")

    resolved, how_counts, missing = [], defaultdict(int), []
    for f in picked["File Name"]:
        path, how = resolve_path(f, by_name, by_uuid)
        how_counts[how] += 1
        if path is None:
            missing.append(f)
        resolved.append(path)

    if how_counts["suffix"] or how_counts["uuid"]:
        print(f"  ファイル名の照合: そのまま一致 {how_counts['name']} 件 / "
              f"接尾辞を除いて一致 {how_counts['suffix']} 件 / "
              f"UUIDで一致 {how_counts['uuid']} 件")
    if missing:
        print(f"  [!] {len(missing)} 件のファイルが見つかりません:")
        for f in missing[:5]:
            print(f"      {f}  (UUID: {file_uuid(f)})")
        sys.exit("ダウンロードが完了しているか、--data-dir が正しいか確認してください。\n"
                 "リネーム済みのファイルはUUID(先頭部分)が残っていれば自動で照合します。")

    picked = picked.assign(_path=resolved)
    n = len(picked)
    sample_ids = picked["Sample ID"].tolist()
    paths = picked["_path"].tolist()
    file_names = [p.name for p in paths]

    first = read_star_file(paths[0])
    genes = first[["gene_id", "gene_name", "gene_type"]].copy()
    ref_ids = first["gene_id"].to_numpy()
    n_genes = len(ref_ids)
    print(f"  基準ファイル: {n_genes:,} 遺伝子 ({file_names[0]})")

    # 事前に配列を確保しておく。merge を繰り返さないので O(n) で済む。
    mats = {v: np.empty((n_genes, n), dtype=np.float64) for v in values}
    for v in values:
        mats[v][:, 0] = first[v].to_numpy(dtype=float)

    started = time.time()
    for i in range(1, n):
        df = read_star_file(paths[i])
        if len(df) != n_genes or not np.array_equal(df["gene_id"].to_numpy(), ref_ids):
            sys.exit(f"遺伝子の並びが基準ファイルと違います: {file_names[i]}\n"
                     f"  (GENCODEのバージョンが混在している可能性があります。"
                     f"1行目の '# gene-model:' を確認してください)")
        for v in values:
            mats[v][:, i] = df[v].to_numpy(dtype=float)
        if (i + 1) % 25 == 0 or i == n - 1:
            el = time.time() - started
            rate = i / el if el else 0
            print(f"    {i + 1:>5} / {n} 検体 ({el:,.0f}秒, 残り約 {(n - i) / rate / 60:,.1f} 分)",
                  flush=True)

    print()
    print("=" * 60)
    print("3. 書き出し")
    print("=" * 60)
    for v in values:
        if args.output and len(values) == 1:
            out = Path(args.output)
        else:
            prefix = args.out_prefix or (Path(args.output).with_suffix("") if args.output
                                         else Path("expression_matrix"))
            out = Path(f"{prefix}_{SUFFIX[v]}.tsv")
        out.parent.mkdir(parents=True, exist_ok=True)
        wide = pd.concat(
            [genes, pd.DataFrame(mats[v], columns=sample_ids)], axis=1
        )
        wide.to_csv(out, sep="\t", index=False)
        print(f"  {v:<20} -> {out}  ({n_genes:,} 遺伝子 × {n} 検体, "
              f"{out.stat().st_size / 1e6:,.0f} MB)")

    print("\npreprocess_tcga.py には次のように渡してください:")
    for v in values:
        flag = {"tpm_unstranded": "--new-tpm", "fpkm_unstranded": "--new-fpkm",
                "fpkm_uq_unstranded": "--new-fpkmuq"}.get(v)
        if flag:
            print(f"    {flag} <上記の {SUFFIX[v]} ファイル>")


if __name__ == "__main__":
    main()
