#!/usr/bin/env python3
"""
make_tnbc_subtype_table.py
==========================
TCGA-BRCA の BCR clinical patient ファイル
(nationwidechildrens_org_clinical_patient_brca.txt)から、免疫組織化学による
ホルモンレセプター/HER2の判定を読み取り、Triple-negative / Non-triple-negative
の2分類CSV(patient_id,subtype)を作る。

出力は preprocess_tcga.py の --tumor-subtype-table にそのまま渡せる。

    python3 make_tnbc_subtype_table.py \
        --input nationwidechildrens_org_clinical_patient_brca.txt \
        --output subtypes/BRCA_TNBC.csv

## 判定ロジック

ER  : er_status_by_ihc
PR  : pr_status_by_ihc
HER2: her2_status_by_ihc を第一とし、それが Positive/Negative 以外
      (Equivocal, Indeterminate, [Not Evaluated], [Not Available])のときだけ
      her2_fish_status で補う。IHCで 2+ = Equivocal ならFISHで確認する、という
      通常の臨床フロー(ASCO/CAP)に合わせている。

3つすべてが Negative      -> Triple-negative
1つでも Positive          -> Non-triple-negative
                             (残りが未確定でも、陽性が1つあれば
                              non-TNであることは確定するため採用する)
上記以外(陽性ゼロだが未確定あり) -> 判定不能。CSVに行を出さない
                             (preprocess_tcga.py 側で自動的に分類対象外になる)

--emit-unknown を付けると、判定不能の患者を 'Unknown' として出力する。
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

POSITIVE = "Positive"
NEGATIVE = "Negative"

ER_COL = "er_status_by_ihc"
PR_COL = "pr_status_by_ihc"
HER2_IHC_COL = "her2_status_by_ihc"
HER2_FISH_COL = "her2_fish_status"


def normalize(value) -> str | None:
    """Positive / Negative のどちらかに確定しているときだけ返す。
    Equivocal, Indeterminate, [Not Evaluated], [Not Available], 空欄 -> None
    """
    v = str(value).strip()
    return v if v in (POSITIVE, NEGATIVE) else None


def read_bcr_clinical(path: str) -> pd.DataFrame:
    """BCRのclinicalファイルはヘッダーが3行ある(1行目=列名、2行目=別名、
    3行目=CDE_ID)。2・3行目を読み飛ばす。"""
    df = pd.read_csv(path, sep="\t", skiprows=[1, 2], low_memory=False)
    if "bcr_patient_barcode" not in df.columns:
        sys.exit(f"{path} に bcr_patient_barcode 列がありません。"
                 f"BCRのclinical patientファイルを指定してください。")
    return df


def classify(row) -> tuple[str | None, str, str, str]:
    er = normalize(row.get(ER_COL))
    pr = normalize(row.get(PR_COL))
    her2 = normalize(row.get(HER2_IHC_COL))
    her2_source = "IHC"
    if her2 is None:
        her2 = normalize(row.get(HER2_FISH_COL))
        her2_source = "FISH" if her2 else "-"

    statuses = [er, pr, her2]
    if all(s == NEGATIVE for s in statuses):
        label = "Triple-negative"
    elif any(s == POSITIVE for s in statuses):
        label = "Non-triple-negative"
    else:
        label = None
    return label, er or "?", pr or "?", f"{her2 or '?'}({her2_source})"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True,
                    help="nationwidechildrens_org_clinical_patient_brca.txt")
    ap.add_argument("--output", default="subtypes/BRCA_TNBC.csv",
                    help="出力CSV(patient_id,subtype)")
    ap.add_argument("--emit-unknown", action="store_true",
                    help="判定不能の患者を 'Unknown' として出力に含める"
                         "(既定は行を出さない=分類対象外)")
    ap.add_argument("--detail-csv",
                    help="患者ごとのER/PR/HER2の内訳も別CSVに書き出す(検証用)")
    args = ap.parse_args()

    df = read_bcr_clinical(args.input)

    for col in (ER_COL, PR_COL, HER2_IHC_COL, HER2_FISH_COL):
        if col not in df.columns:
            sys.exit(f"必要な列 '{col}' が見つかりません。列名: {list(df.columns)[:20]} ...")

    records, details = [], []
    for _, row in df.iterrows():
        pid = str(row["bcr_patient_barcode"]).strip()
        label, er, pr, her2 = classify(row)
        details.append({"patient_id": pid, "ER": er, "PR": pr, "HER2": her2,
                        "subtype": label or "Unknown"})
        if label is None and not args.emit_unknown:
            continue
        records.append({"patient_id": pid, "subtype": label or "Unknown"})

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records).to_csv(out, index=False)

    detail_df = pd.DataFrame(details)
    counts = detail_df["subtype"].value_counts()
    print(f"読み込んだ患者: {len(df)} 名")
    for name, n in counts.items():
        note = ""
        if name == "Unknown":
            note = " (陽性ゼロだが未確定の項目あり" + \
                   ("。--emit-unknown により出力に含めています)" if args.emit_unknown
                    else "。出力には含めていません)")
        print(f"  {name}: {n}{note}")

    n_fish = detail_df["HER2"].str.contains(r"\(FISH\)").sum()
    print(f"HER2をFISHで補った患者: {n_fish} 名")
    print(f"出力: {out} ({len(records)} 行)")

    if args.detail_csv:
        Path(args.detail_csv).parent.mkdir(parents=True, exist_ok=True)
        detail_df.to_csv(args.detail_csv, index=False)
        print(f"内訳: {args.detail_csv}")


if __name__ == "__main__":
    main()
