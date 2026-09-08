#!/usr/bin/env python3
"""
make_clinical_subtype_table.py
==============================
BCRの臨床ファイル(または論文の補足表)の任意の列から、
preprocess_tcga.py に渡せる patient_id,subtype 形式のCSVを作る。

組織型のように「元から表に入っている分類」は、論文の補足資料を探すより
臨床ファイルを見た方が早いことが多い。

## まず中身を調べる

    python3 scripts/make_clinical_subtype_table.py \\
        --input nationwidechildrens_org_clinical_patient_stad.txt --inspect

列名と値の分布が出るので、使えそうな列を選ぶ。

## 変換する

値をそのまま使う場合:

    python3 scripts/make_clinical_subtype_table.py \\
        --input nationwidechildrens_org_clinical_patient_stad.txt \\
        --column histological_type \\
        --output subtypes/STAD_histology.csv

値をまとめ直す場合(Lauren分類など):

    python3 scripts/make_clinical_subtype_table.py \\
        --input nationwidechildrens_org_clinical_patient_stad.txt \\
        --column histological_type \\
        --preset lauren \\
        --output subtypes/STAD_Lauren.csv

--preset lauren は、値に "diffuse" / "intestinal" / "signet ring" などが
含まれるかを見て Diffuse / Intestinal / Mixed に振り分ける。
実際の表記は必ず --inspect で確認し、想定と違えば --map で明示指定すること。

任意の対応を指定する場合:

    --map "Stomach, Adenocarcinoma, Diffuse Type=Diffuse" \\
    --map "Stomach, Intestinal Adenocarcinoma, Tubular Type=Intestinal"

--map を指定すると、対応表に無い値は出力されない(分類対象外になる)。
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

BARCODE_RE = re.compile(r"(TCGA-\w{2}-\w{4})")
BARCODE_PAT = r"TCGA-\w{2}-\w{4}"
BARCODE_CANDIDATES = [
    "bcr_patient_barcode", "patient_id", "case id", "sample id",
    "tumor sample id", "submitter_id", "patient",
]
NULLS = {"", "na", "n/a", "nan", "none", "-", "--", "null",
         "[not available]", "[not applicable]", "[unknown]",
         "[discrepancy]", "not available", "unknown", "#n/a"}

# 組織型の記載からLauren分類を推定するための規則。
# 判定は「上から順に最初に当たったもの」を採用する(印環細胞癌はびまん型に
# 含めるのが一般的なので、diffuse より先に置いてある)。
LAUREN_RULES = [
    ("Diffuse",    ["signet ring", "diffuse", "poorly cohesive"]),
    ("Intestinal", ["intestinal", "tubular", "papillary"]),
    ("Mixed",      ["mixed"]),
]

PRESETS = {"lauren": LAUREN_RULES}


def norm(v) -> str | None:
    s = str(v).strip()
    return None if s.lower() in NULLS else s


def read_table(path: str) -> pd.DataFrame:
    """BCRの臨床ファイルは列名の下に説明行が1〜2行入る。
    TCGAバーコードが最初に現れる行を「最初のデータ行」とみなし、
    その1つ上をヘッダーとして読み直す。"""
    sep = "\t" if Path(path).suffix.lower() in (".txt", ".tsv") else ","
    if Path(path).suffix.lower() in (".xlsx", ".xls"):
        raw = pd.read_excel(path, header=None, dtype=str)
        reader = lambda hdr: pd.read_excel(path, header=hdr, dtype=str)
    else:
        raw = pd.read_csv(path, sep=sep, header=None, dtype=str,
                          low_memory=False, nrows=30)
        reader = lambda hdr: pd.read_csv(path, sep=sep, header=hdr, dtype=str,
                                         low_memory=False)

    first_data = None
    for i in range(min(30, len(raw))):
        if raw.iloc[i].map(lambda v: bool(BARCODE_RE.search(str(v)))).any():
            first_data = i
            break
    if first_data is None:
        sys.exit(f"{path} にTCGAバーコードらしき値が見つかりません。")

    # 説明行の入り方が2通りあるので区別する。
    #   BCRの臨床ファイル : 1行目が列名、その下に別名行・CDE_ID行が続く
    #                       -> ヘッダーは0行目、間の行を読み飛ばす
    #   論文の補足表      : 表題や説明が先にあり、その下に列名の行が来る
    #                       -> ヘッダーは「最初のデータ行の1つ上」
    has_cde = any("CDE_ID" in " ".join(str(v) for v in raw.iloc[i].tolist())
                  for i in range(min(first_data, len(raw))))
    if has_cde or first_data <= 1:
        header_row, skip = 0, list(range(1, first_data))
    else:
        header_row, skip = first_data - 1, []

    df = reader(header_row)
    if skip:
        # ヘッダー直下の説明行を落とす(読み込み後に位置で除く)
        df = df.drop(index=[i - 1 for i in skip if i - 1 < len(df)])
    if header_row or skip:
        n_skipped = header_row + len(skip)
        print(f"  (説明行 {n_skipped} 行を読み飛ばしました)")
    # ヘッダーの直下に残った説明行(バーコードを含まない行)を落とす
    bc_col = find_barcode_column(df)
    if bc_col is not None:
        df = df[df[bc_col].astype(str).str.contains(BARCODE_PAT, na=False, regex=True)]
    return df.reset_index(drop=True)


def find_barcode_column(df: pd.DataFrame):
    lower = {str(c).strip().lower(): c for c in df.columns}
    for cand in BARCODE_CANDIDATES:
        if cand in lower:
            return lower[cand]
    best, best_n = None, 0
    for c in df.columns:
        n = df[c].astype(str).str.contains(BARCODE_PAT, na=False, regex=True).sum()
        if n > best_n:
            best, best_n = c, n
    return best if best_n >= 5 else None


def inspect(df: pd.DataFrame) -> None:
    print(f"\n{len(df)} 行 × {len(df.columns)} 列")
    print(f"患者ID列の推定: {find_barcode_column(df)!r}\n")
    print("分類に使えそうな列(値が2〜20種類のもの):")
    shown = 0
    for c in df.columns:
        vals = df[c].map(norm).dropna()
        uniq = vals.unique()
        if 2 <= len(uniq) <= 20:
            counts = vals.value_counts()
            detail = ", ".join(f"{k}({v})" for k, v in list(counts.items())[:8])
            if len(counts) > 8:
                detail += f", ... 他{len(counts)-8}種"
            print(f"  {str(c)[:40]:<42}{detail}")
            shown += 1
    if shown == 0:
        print("  (該当なし。--inspect-all で全列を表示できます)")


def apply_rules(value: str, rules) -> str | None:
    low = value.lower()
    for label, keys in rules:
        if any(k in low for k in keys):
            return label
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True)
    ap.add_argument("--inspect", action="store_true",
                    help="列名と値の分布を表示して終了する")
    ap.add_argument("--inspect-all", action="store_true",
                    help="値の種類が多い列も含めてすべて表示する")
    ap.add_argument("--column", help="分類に使う列名")
    ap.add_argument("--barcode-column", help="患者ID列を明示指定する")
    ap.add_argument("--output", help="出力CSV(patient_id,subtype)")
    ap.add_argument("--preset", choices=sorted(PRESETS),
                    help="値のまとめ方の既定規則")
    ap.add_argument("--map", action="append", default=[], metavar="元の値=表示名",
                    help="値の対応を明示指定する。複数回指定可。"
                         "指定した場合、対応表に無い値は出力しない")
    ap.add_argument("--min-n", type=int, default=1,
                    help="この人数に満たない群を出力から外す(既定 1 = 外さない)")
    args = ap.parse_args()

    df = read_table(args.input)

    if args.inspect or args.inspect_all:
        if args.inspect_all:
            for c in df.columns:
                vals = df[c].map(norm).dropna()
                print(f"  {str(c)[:40]:<42}[{vals.nunique()} 種類]")
        else:
            inspect(df)
        print("\n次は --column で列を選んで変換してください。")
        return

    if not args.column or not args.output:
        sys.exit("--column と --output を指定してください(--inspect で列を確認できます)")
    if args.column not in df.columns:
        sys.exit(f"列 '{args.column}' がありません。--inspect で確認してください。")

    bc = args.barcode_column or find_barcode_column(df)
    if not bc:
        sys.exit("患者ID列を特定できませんでした。--barcode-column で指定してください。")
    print(f"患者ID列: {bc!r} / 分類に使う列: {args.column!r}")

    explicit = {}
    for spec in args.map:
        if "=" not in spec:
            sys.exit(f"--map は '元の値=表示名' の形式です: {spec}")
        k, v = spec.split("=", 1)
        explicit[k.strip().lower()] = v.strip()
    rules = PRESETS.get(args.preset)

    records, unmatched = [], {}
    for _, r in df.iterrows():
        m = BARCODE_RE.search(str(r[bc]))
        value = norm(r[args.column])
        if not m or not value:
            continue
        if explicit:
            label = explicit.get(value.lower())
        elif rules:
            label = apply_rules(value, rules)
        else:
            label = value
        if label is None:
            unmatched[value] = unmatched.get(value, 0) + 1
            continue
        records.append({"patient_id": m.group(1), "subtype": label})

    out = pd.DataFrame(records, columns=["patient_id", "subtype"])
    if args.min_n > 1:
        counts = out["subtype"].value_counts()
        keep = set(counts[counts >= args.min_n].index)
        dropped = set(counts.index) - keep
        if dropped:
            print(f"  {args.min_n} 人未満の群を除外: {sorted(dropped)}")
        out = out[out["subtype"].isin(keep)]

    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, index=False)

    print(f"\n出力: {path} ({len(out)} 患者)")
    for k, v in out["subtype"].value_counts().items():
        print(f"    {k}: {v}")
    if unmatched:
        print("\n  [注意] 対応が付かず出力しなかった値:")
        for k, v in sorted(unmatched.items(), key=lambda x: -x[1]):
            print(f"      {k!r}: {v} 件")
        print("      必要なら --map で明示指定してください。")


if __name__ == "__main__":
    main()
