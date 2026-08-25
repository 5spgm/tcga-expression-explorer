#!/usr/bin/env python3
"""
make_icluster_subtype_table.py
===============================
Hoadley et al. 2018 (Cell 173(2):291-304.e6) の Table S6 (mmc6.xlsx, iCluster
membership) を、preprocess_tcga.py が読める汎用フォーマット
(patient_id, subtype の2列CSV)に変換する。

Table S6 の中身は「Sample ID (12桁の患者バーコード), iCluster (1-28の番号)」
の2列のみで、論文本文にある `C4:pan-GI (CRC)` のような名前は別途対応させる
必要がある。本スクリプトはその対応表 (ICLUSTER_LABELS) を内蔵している。

## 使い方
    python3 make_icluster_subtype_table.py \\
        --input mmc6.xlsx \\
        --output icluster_subtypes.csv

出力CSVは patient_id,subtype の2列。この後 preprocess_tcga.py に
--tumor-subtype-table icluster_subtypes.csv として渡す。

## 注意
論文本文中で名前が言及されていた iCluster (C1-C16, C18-C21, C23-C28) には
記述的な名前を付けている。C17, C22 は本文中に名前の記載が見当たらなかった
ため、"iCluster 17" のような番号のみの表記になる。正確な名前が必要な場合は
論文 Figure 2 や Table S6 の元シート(色分け凡例)を直接確認すること。
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

# Hoadley et al. 2018 本文 (Results / Discussion) に記載されている名称。
# 該当箇所が見つからなかった C17, C22 は番号のみのプレースホルダ。
ICLUSTER_LABELS = {
    1: "C1:STAD (EBV-CIMP)",
    2: "C2:BRCA (HER2 amp)",
    3: "C3:mesenchymal (immune)",
    4: "C4:pan-GI (CRC)",
    5: "C5:CNS/endocrine",
    6: "C6:OV",
    7: "C7:mixed (Chr9 del)",
    8: "C8:UCEC",
    9: "C9:ACC/KICH",
    10: "C10:pan-SCC",
    11: "C11:LGG (IDH1 mut)",
    12: "C12:THCA",
    13: "C13:mixed (Chr8 del)",
    14: "C14:LUAD",
    15: "C15:SKCM/UVM",
    16: "C16:PRAD",
    17: "C17",
    18: "C18:pan-GI (MSI)",
    19: "C19:BRCA (luminal)",
    20: "C20:mixed (stromal/immune)",
    21: "C21:DLBC",
    22: "C22",
    23: "C23:GBM/LGG (IDH1wt)",
    24: "C24:LAML",
    25: "C25:pan-SCC (Chr11 amp)",
    26: "C26:LIHC",
    27: "C27:pan-SCC (HPV)",
    28: "C28:pan-kidney",
}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="mmc6.xlsx (Table S6) のパス")
    ap.add_argument("--output", required=True, help="出力CSVのパス")
    args = ap.parse_args()

    df = pd.read_excel(args.input, header=1)
    if df.shape[1] < 2:
        sys.exit(f"想定外の列数です: {df.shape}")

    df = df.iloc[:, :2]
    df.columns = ["patient_id", "icluster"]
    df = df.dropna(subset=["patient_id", "icluster"])
    df["icluster"] = df["icluster"].astype(int)
    df["subtype"] = df["icluster"].map(ICLUSTER_LABELS).fillna(
        df["icluster"].apply(lambda n: f"C{n}")
    )

    out = df[["patient_id", "subtype"]]
    out.to_csv(args.output, index=False)
    print(f"{len(out)} 件の患者IDを {args.output} に出力しました")
    print(out["subtype"].value_counts())


if __name__ == "__main__":
    main()
