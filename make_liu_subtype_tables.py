#!/usr/bin/env python3
"""
make_liu_subtype_tables.py
===========================
Liu et al. 2018 (Cell 173(4):963-985.e16, "Comparative Molecular Analysis of
Gastrointestinal Adenocarcinomas") の Table S1 (mmc2.xlsx, Master Patient
Table) から、複数のsubtype分類をそれぞれ patient_id,subtype の2列CSVとして
書き出す。

出力される分類(既定):
  - Molecular_Subtype        論文のメイン分類 (CIN / MSI / GS / HM-SNV)
  - MSI_Status                マイクロサテライト不安定性 (MSS / MSI-H / MSI-L)
  - CIMP                      Hypermethylation category (CIMP-H / CRC CIMP-L / Non-CIMP / GEA CIMP-L)
  - Colorectal_CMS            大腸がんのCMS分類 (Guinney et al., 2015; COAD/READのみ)

## 使い方
    python3 make_liu_subtype_tables.py \\
        --input mmc2.xlsx \\
        --out-dir ./subtypes

    -> ./subtypes/Molecular_Subtype.csv, MSI_Status.csv, CIMP.csv, Colorectal_CMS.csv
       が生成される。それぞれ preprocess_tcga.py の --subtype-scheme に渡す。
"""

import argparse
from pathlib import Path

import pandas as pd

# 出力ファイル名 -> (元の列名, 表示ラベルの接頭辞)
COLUMN_MAP = {
    "Molecular_Subtype": "Molecular_Subtype",
    "MSI_Status": "MSI Status",
    "CIMP": "Hypermethylation category",
    "Colorectal_CMS": "Colorectal CMS",
}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="mmc2.xlsx (Table S1) のパス")
    ap.add_argument("--out-dir", required=True, help="CSV出力先ディレクトリ")
    ap.add_argument("--sheet", default="Master Patient Table")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_excel(args.input, sheet_name=args.sheet, header=1)
    if "TCGA Participant Barcode" not in df.columns:
        raise SystemExit(f"'TCGA Participant Barcode' 列が見つかりません。列一覧: {list(df.columns)[:10]}...")

    for out_name, src_col in COLUMN_MAP.items():
        if src_col not in df.columns:
            print(f"  [skip] 列 '{src_col}' が見つからないためスキップ: {out_name}")
            continue
        sub = df[["TCGA Participant Barcode", src_col]].dropna()
        sub.columns = ["patient_id", "subtype"]
        sub = sub[sub["subtype"].astype(str).str.strip() != ""]
        out_path = out_dir / f"{out_name}.csv"
        sub.to_csv(out_path, index=False)
        print(f"[{out_name}] {len(sub)} 件 -> {out_path}")
        print(f"  内訳: {sub['subtype'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
