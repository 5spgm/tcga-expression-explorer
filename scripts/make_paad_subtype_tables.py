#!/usr/bin/env python3
"""
make_paad_subtype_tables.py
===========================
TCGA-PAAD marker paper (Raphael et al., Cancer Cell 2017) の補足表から、
preprocess_tcga.py に渡せる patient_id,subtype 形式のCSVを作る。

膵臓がん特有の分類:
  Moffitt   : Basal-like / Classical          (予後との関連が最も確立)
  Collisson : Classical / Exocrine-like / QM-PDA
  Bailey    : Squamous / Immunogenic / Progenitor / ADEX

加えて、**PDACでない検体の除外リスト**も作れる。TCGA-PAADには神経内分泌腫瘍・
腺房細胞癌・IPMN由来などが混入しており、論文でも解析から外されている。
膵管腺癌だけを見たい場合は、生成した除外リストを preprocess_tcga.py の
--exclude-samples に渡すこと。

## まず中身を調べる

補足表のレイアウトは論文・版によって違うため、最初に --inspect で
シート名・列名・値の分布を確認する:

    python3 make_paad_subtype_tables.py --input mmc2.xlsx --inspect

## 変換する

自動検出に任せる場合:

    python3 make_paad_subtype_tables.py \\
        --input mmc2.xlsx --out-dir ./subtypes/paad

列名が自動検出できない場合は明示指定する:

    python3 make_paad_subtype_tables.py \\
        --input mmc2.xlsx --out-dir ./subtypes/paad \\
        --sheet "Table S1" \\
        --barcode-column "Tumor Sample ID" \\
        --scheme "Moffitt=mRNA Moffitt clusters" \\
        --scheme "Bailey=mRNA Bailey Clusters"
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

BARCODE_CANDIDATES = [
    "tumor sample id", "patient id", "bcr_patient_barcode", "case id",
    "sample id", "tcga id", "submitter_id", "patient", "sampleid",
]

# 列名の一部にこれが含まれていればその分類とみなす(小文字で比較)。
# 値はそのまま出さず、下の CODE_LABELS で意味のある名前に置き換える。
SCHEME_PATTERNS = {
    "Moffitt": ["moffitt"],
    "Collisson": ["collisson"],
    "Bailey": ["bailey"],
    # 純度そのものを分類方式として出す。膵臓がんは間質が主体になりやすく、
    # 解析対象150検体でも high 76 / low 74 と半々。Moffittなどの発現ベースの
    # 分類が純度に引きずられていないかを、その場で確認できるようにする。
    "Purity Class": ["purity class"],
    "Copy Number": ["copy number cluster"],
    "Methylation": ["methylation cluster"],
    "miRNA": ["mirna cluster"],
    "lncRNA": ["lncrna cluster"],
    "RPPA": ["rppa cluster"],
}

# Cancer Cell 2017 の Table S1 では、分類が 1/2/3/4 の数値コードで入っており、
# 意味は列名に埋め込まれている
#   例: 'mRNA Moffitt clusters (All 150 Samples) 1basal  2classical'
# そのまま出すとサイト上で "1" "2" としか表示されず解読できないので、
# ここで意味のあるラベルに置き換える。
CODE_LABELS = {
    "Moffitt":   {"1": "Basal-like", "2": "Classical"},
    "Collisson": {"1": "Classical", "2": "Exocrine-like", "3": "QM-PDA"},
    "Bailey":    {"1": "Squamous", "2": "Immunogenic",
                  "3": "Progenitor", "4": "ADEX"},
    # high/low だけだと何の高低か分からないので明示する
    "Purity Class": {"high": "High purity", "low": "Low purity"},
}

# 分類として意味を持たない値(低純度のため判定不能、など)
NON_LABELS = {"lowpurity", "low purity", "not applicable", "nc", "unclassified"}

# 「全150検体版」と「高純度76検体のみ版」が両方ある場合、既定では
# 前者を採る(検体数が多く、サイト上で群として成立するため)。
PREFER_ALL_SAMPLES = "all 150 samples"
AVOID_SUBSET = "high purity samples only"

# 非PDACの判定に使いそうな列
HISTOLOGY_PATTERNS = ["histolog", "diagnosis", "pathology", "tumor type", "included"]

BARCODE_RE = re.compile(r"(TCGA-\w{2}-\w{4})")
# str.contains 用(キャプチャグループなし。あると警告が出る)
BARCODE_PAT = r"TCGA-\w{2}-\w{4}"
NULLS = {"", "na", "n/a", "nan", "none", "-", "--", "null", "not available",
         "[not available]", "[not applicable]", "not applicable", "#n/a"}


def norm(v) -> str | None:
    s = str(v).strip()
    return None if s.lower() in NULLS else s


def read_sheet(path: str, sheet: str | None) -> tuple[pd.DataFrame, str]:
    """ヘッダー行がファイルの先頭にあるとは限らない(説明行が数行入ることがある)。
    TCGAバーコードらしき値が最も多く現れる行をヘッダーとみなして読み直す。"""
    xl = pd.ExcelFile(path)
    name = sheet or xl.sheet_names[0]
    if name not in xl.sheet_names:
        sys.exit(f"シート '{name}' がありません。存在するシート: {xl.sheet_names}")

    raw = pd.read_excel(xl, sheet_name=name, header=None, dtype=str)
    # TCGAバーコードが最初に現れる行が「最初のデータ行」。その1つ上がヘッダー。
    first_data_row = None
    for i in range(min(30, len(raw))):
        if raw.iloc[i].map(lambda v: bool(BARCODE_RE.search(str(v)))).any():
            first_data_row = i
            break
    if first_data_row is None:
        sys.exit(f"シート '{name}' にTCGAバーコードらしき値が見つかりません。"
                 f"--sheet で別のシートを指定してください。"
                 f"(存在するシート: {xl.sheet_names})")
    header_row = max(0, first_data_row - 1)
    df = pd.read_excel(xl, sheet_name=name, header=header_row, dtype=str)
    df = df.dropna(how="all").dropna(axis=1, how="all")
    if header_row:
        print(f"  (先頭 {header_row} 行は説明行と判断して読み飛ばしました)")
    return df, name


def find_barcode_column(df: pd.DataFrame) -> str | None:
    lower = {str(c).strip().lower(): c for c in df.columns}
    for cand in BARCODE_CANDIDATES:
        if cand in lower:
            return lower[cand]
    # 名前で見つからなければ、TCGAバーコードが最も多い列を選ぶ
    best, best_n = None, 0
    for c in df.columns:
        n = df[c].astype(str).str.contains(BARCODE_PAT, na=False, regex=True).sum()
        if n > best_n:
            best, best_n = c, n
    return best if best_n >= 5 else None


def inspect(df: pd.DataFrame, sheet_name: str) -> None:
    print(f"\nシート '{sheet_name}': {len(df)} 行 × {len(df.columns)} 列\n")
    bc = find_barcode_column(df)
    print(f"バーコード列の推定: {bc!r}\n")
    print("列の一覧(カテゴリ的な列は値の分布も表示):")
    for c in df.columns:
        vals = df[c].map(norm).dropna()
        uniq = vals.unique()
        marker = ""
        low = str(c).strip().lower()
        for scheme, pats in SCHEME_PATTERNS.items():
            if any(p in low for p in pats):
                marker = f"  ← {scheme} らしい"
        if any(p in low for p in HISTOLOGY_PATTERNS):
            marker = "  ← 組織型/除外判定に使えそう"
        if 1 < len(uniq) <= 12:
            counts = vals.value_counts()
            detail = ", ".join(f"{k}({v})" for k, v in counts.items())
            print(f"  {str(c)[:46]:<48}{detail}{marker}")
        else:
            print(f"  {str(c)[:46]:<48}[{len(uniq)} 種類の値]{marker}")


def write_csv(records, path: Path) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records, columns=["patient_id", "subtype"]).to_csv(path, index=False)
    return len(records)


def write_exclude_lists_from_sheets(path: str, out_dir: Path) -> None:
    """シート名から除外系のシートを探し、患者IDの一覧を書き出す。"""
    xl = pd.ExcelFile(path)
    targets = {
        "PAAD_exclude_non_pdac.txt": (["exclud"], "膵管腺癌でないと判定された検体"),
        "PAAD_exclude_pseudonormal.txt": (["pseudonormal", "pseudo normal", "low purity"],
                                          "腫瘍細胞含有率が極端に低い検体"),
    }
    produced = {}
    for fname, (pats, desc) in targets.items():
        sheets = [s for s in xl.sheet_names
                  if any(p in s.lower() for p in pats)]
        if not sheets:
            continue
        ids = set()
        for sheet in sheets:
            d = pd.read_excel(xl, sheet_name=sheet, header=1, dtype=str)
            for col in d.columns:
                for v in d[col].dropna().astype(str):
                    m = BARCODE_RE.search(v)
                    if m:
                        ids.add(m.group(1))
                if ids:
                    break   # 最初にバーコードが見つかった列だけ使う
        if not ids:
            continue
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / fname).write_text("\n".join(sorted(ids)) + "\n")
        produced[fname] = (len(ids), desc, sheets)

    if not produced:
        return
    print("\n除外リスト(シート由来):")
    all_ids = set()
    for fname, (n, desc, sheets) in produced.items():
        print(f"  {desc}: {n} 患者 -> {out_dir / fname}")
        print(f"      (シート: {', '.join(sheets)})")
        all_ids |= set((out_dir / fname).read_text().split())
    if len(produced) > 1:
        (out_dir / "PAAD_exclude_all.txt").write_text("\n".join(sorted(all_ids)) + "\n")
        print(f"  上記をまとめたもの: {len(all_ids)} 患者 -> "
              f"{out_dir / 'PAAD_exclude_all.txt'}")
    print("  preprocess_tcga.py の --exclude-samples に渡してください")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="補足表 (xlsx)")
    ap.add_argument("--sheet", help="シート名(省略時は先頭シート)")
    ap.add_argument("--inspect", action="store_true",
                    help="シート名・列名・値の分布を表示して終了する")
    ap.add_argument("--out-dir", default="./subtypes/paad")
    ap.add_argument("--barcode-column", help="患者IDの列を明示指定する")
    ap.add_argument("--scheme", action="append", default=[], metavar="表示名=列名",
                    help="分類を明示指定する。複数回指定可")
    ap.add_argument("--exclude-column",
                    help="非PDAC検体の判定に使う列(組織型など)")
    ap.add_argument("--exclude-values",
                    help="--exclude-column のうち除外する値をカンマ区切りで指定")
    ap.add_argument("--no-exclude-lists", action="store_true",
                    help="ExcludedSamples / PseudoNormals シートからの"
                         "除外リスト生成を行わない")
    args = ap.parse_args()

    df, sheet_name = read_sheet(args.input, args.sheet)

    if args.inspect:
        inspect(df, sheet_name)
        print("\n次は、上の一覧を見て --scheme '表示名=列名' で変換してください。")
        return

    bc = args.barcode_column or find_barcode_column(df)
    if not bc:
        sys.exit("患者ID列を特定できませんでした。--inspect で確認し、"
                 "--barcode-column で指定してください。")
    print(f"患者ID列: {bc!r}")

    # 患者IDは TCGA-XX-XXXX に切り詰める(検体バーコードで入っている場合に対応)
    df["_pid"] = df[bc].astype(str).map(
        lambda v: (BARCODE_RE.search(v).group(1) if BARCODE_RE.search(v) else None))
    df = df[df["_pid"].notna()]
    print(f"患者ID を持つ行: {len(df)}")

    # --- 分類を決める --------------------------------------------------
    schemes: dict[str, str] = {}
    for spec in args.scheme:
        if "=" not in spec:
            sys.exit(f"--scheme は '表示名=列名' の形式で指定してください: {spec}")
        name, col = spec.split("=", 1)
        if col.strip() not in df.columns:
            sys.exit(f"列 '{col.strip()}' がありません。--inspect で確認してください。")
        schemes[name.strip()] = col.strip()

    if not schemes:  # 自動検出
        # 同じ分類に複数の列(全検体版 / 高純度のみ版)がある場合は全検体版を優先
        for prefer in (True, False):
            for c in df.columns:
                low = str(c).strip().lower()
                is_all = PREFER_ALL_SAMPLES in low
                if prefer != is_all:
                    continue
                for scheme, pats in SCHEME_PATTERNS.items():
                    if any(p in low for p in pats) and scheme not in schemes:
                        if not prefer and AVOID_SUBSET in low:
                            continue  # 高純度のみ版は、全検体版が無いときだけ使う
                        schemes[scheme] = c
        if schemes:
            print("\n自動検出した分類:")
            for k, v in schemes.items():
                print(f"  {k:<22} <- 列 {v!r}")
        else:
            sys.exit("分類の列を自動検出できませんでした。--inspect で確認し、"
                     "--scheme で指定してください。")

    out_dir = Path(args.out_dir)
    print()
    for name, col in schemes.items():
        code_map = CODE_LABELS.get(name, {})
        records, n_dropped = [], 0
        for _, r in df.iterrows():
            label = norm(r[col])
            if not label:
                continue
            if label.strip().lower() in NON_LABELS:
                n_dropped += 1        # 低純度で判定不能。分類対象外にする
                continue
            mapped = code_map.get(label)
            if mapped is None:
                # 意味の分かっていない数値クラスタは "Cluster 1" の形にする
                # (サイト上で "1" とだけ出ると何のことか分からないため)
                mapped = f"Cluster {label}" if label.isdigit() else label
            records.append({"patient_id": r["_pid"], "subtype": mapped})
        if n_dropped:
            print(f"  ({name}: 判定不能 {n_dropped} 件は出力しません)")
        if not records:
            print(f"  [スキップ] {name}: 有効な値がありませんでした")
            continue
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", name)
        path = out_dir / f"PAAD_{safe}.csv"
        write_csv(records, path)
        counts = pd.Series([r["subtype"] for r in records]).value_counts()
        print(f"  {name}: {len(records)} 患者 -> {path}")
        print(f"      {", ".join(f"{k}: {v}" for k, v in counts.items())}")

    # --- 別シートからの除外リスト --------------------------------------
    # Cancer Cell 2017 の Table S1 には、解析から外された検体が別シートに
    # まとめられている:
    #   ExcludedSamples            膵管腺癌ではないと判定されたもの
    #   PseudoNormals (low purity) 腫瘍細胞含有率が極端に低いもの
    # 膵臓がんは間質が主体になりやすく、後者は「腫瘍」として扱うと
    # 実質的に間質の発現を見ていることになるため、用途によっては外す。
    if not args.no_exclude_lists:
        write_exclude_lists_from_sheets(args.input, out_dir)

    if args.exclude_column:
        if args.exclude_column not in df.columns:
            sys.exit(f"列 '{args.exclude_column}' がありません。")
        drop_vals = {v.strip().lower()
                     for v in (args.exclude_values or "").split(",") if v.strip()}
        if not drop_vals:
            sys.exit("--exclude-values も指定してください "
                     "(--inspect で値の一覧を確認できます)")
        hit = df[df[args.exclude_column].astype(str).str.strip().str.lower().isin(drop_vals)]
        path = out_dir / "PAAD_exclude_non_pdac.txt"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(sorted(set(hit["_pid"]))) + "\n")
        print(f"\n  非PDAC等の除外対象: {hit['_pid'].nunique()} 患者 -> {path}")
        print(f"      preprocess_tcga.py の --exclude-samples に渡してください")


if __name__ == "__main__":
    main()
