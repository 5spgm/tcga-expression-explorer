#!/usr/bin/env python3
"""
make_immune_subtype_tables.py
=============================
Thorsson et al., Immunity 2018 "The Immune Landscape of Cancer" の
Supplementary Table (mmc2.xlsx, シート 'PanImmune_MS') から、
preprocess_tcga.py に渡せる patient_id,subtype 形式のCSVを作る。

このファイル1つから**2種類**の分類が取り出せる。

1. Immune Subtype (C1-C6) — 全33がん種横断
   'Immune Subtype' 列。Hoadleyのi Clusterと同じく汎がん種なので、
   どのがん種でもそのまま使い回せる。
   C1 Wound Healing / C2 IFN-gamma Dominant / C3 Inflammatory /
   C4 Lymphocyte Depleted / C5 Immunologically Quiet / C6 TGF-beta Dominant

2. TCGA Subtype — がん種ごとの既発表サブタイプ
   'TCGA Subtype' 列。`BRCA.LumA` のように「がん種.ラベル」形式で入っている。
   **乳がんではこれがPAM50** (LumA / LumB / Basal / Her2 / Normal)。
   接頭辞ごとに別ファイルへ切り出す。

使い方:

    python3 make_immune_subtype_tables.py \
        --input 1-s2_0-S1074761318301213-mmc2.xlsx \
        --out-dir ./subtypes

生成物:
    subtypes/Immune_Subtype.csv           全がん種共通(そのまま使い回せる)
    subtypes/tcga_subtype/BRCA.csv        PAM50
    subtypes/tcga_subtype/GI.csv          消化管腺癌(COAD/READ/STAD/ESCAを含む)
    subtypes/tcga_subtype/<接頭辞>.csv    …各がん種

注意: 'TCGA Subtype' の接頭辞は必ずしもTCGAのstudy略号と一致しない
(`GBM_LGG`, `GI`, `OVCA` など複数がん種をまとめた単位がある)。
どのがん種にどのファイルを渡すかは、生成時に表示される対応表で確認すること。
該当患者がいないがん種に渡しても、単にプルダウンに出てこないだけで害はない。
"""

import argparse
import re
from pathlib import Path

import pandas as pd

SHEET = "PanImmune_MS"
BARCODE_COL = "TCGA Participant Barcode"
STUDY_COL = "TCGA Study"
IMMUNE_COL = "Immune Subtype"
TCGA_SUBTYPE_COL = "TCGA Subtype"

# 論文中の呼称。C1..C6 だけでは中身が分からないので、表示名に併記する。
IMMUNE_NAMES = {
    "C1": "C1 Wound Healing",
    "C2": "C2 IFN-gamma Dominant",
    "C3": "C3 Inflammatory",
    "C4": "C4 Lymphocyte Depleted",
    "C5": "C5 Immunologically Quiet",
    "C6": "C6 TGF-beta Dominant",
}

SAFE = re.compile(r"[^A-Za-z0-9_.-]")


def write_csv(records, path: Path) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records, columns=["patient_id", "subtype"]).to_csv(path, index=False)
    return len(records)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="mmc2.xlsx (Thorsson 2018)")
    ap.add_argument("--out-dir", default="./subtypes", help="出力先ディレクトリ")
    ap.add_argument("--raw-immune-labels", action="store_true",
                    help="Immune Subtypeを 'C1' のような素のラベルで出力する"
                         "(既定は 'C1 Wound Healing' のように説明を併記)")
    args = ap.parse_args()

    df = pd.read_excel(args.input, sheet_name=SHEET)
    for col in (BARCODE_COL, STUDY_COL, IMMUNE_COL, TCGA_SUBTYPE_COL):
        if col not in df.columns:
            raise SystemExit(f"列 '{col}' がありません。列名: {list(df.columns)[:10]} ...")

    out_dir = Path(args.out_dir)
    print(f"読み込み: {len(df)} 行 / {df[STUDY_COL].nunique()} がん種\n")

    # ---- 1. Immune Subtype (全がん種共通) ----------------------------
    imm = df[[BARCODE_COL, IMMUNE_COL]].dropna()
    records = [
        {"patient_id": str(r[BARCODE_COL]).strip(),
         "subtype": str(r[IMMUNE_COL]).strip() if args.raw_immune_labels
                    else IMMUNE_NAMES.get(str(r[IMMUNE_COL]).strip(), str(r[IMMUNE_COL]).strip())}
        for _, r in imm.iterrows()
    ]
    path = out_dir / "Immune_Subtype.csv"
    n = write_csv(records, path)
    print(f"[1] Immune Subtype -> {path} ({n} 患者)")
    counts = imm[IMMUNE_COL].value_counts().sort_index()
    for k, v in counts.items():
        print(f"      {IMMUNE_NAMES.get(k, k)}: {v}")

    # ---- 2. TCGA Subtype (接頭辞ごとに分割) ---------------------------
    print(f"\n[2] TCGA Subtype -> {out_dir / 'tcga_subtype'}/")
    sub = df[[BARCODE_COL, STUDY_COL, TCGA_SUBTYPE_COL]].dropna(subset=[TCGA_SUBTYPE_COL]).copy()
    parts = sub[TCGA_SUBTYPE_COL].astype(str).str.split(".", n=1, expand=True)
    sub["_prefix"] = parts[0].str.strip()
    sub["_label"] = (parts[1] if parts.shape[1] > 1 else "").fillna("").str.strip()
    # ラベルが空('GBM_LGG.' など)や '-' の行は分類不能なので落とす
    sub = sub[~sub["_label"].isin(["", "-", "nan"])]

    summary = []
    for prefix, grp in sub.groupby("_prefix"):
        records = [{"patient_id": str(r[BARCODE_COL]).strip(), "subtype": r["_label"]}
                   for _, r in grp.iterrows()]
        fname = SAFE.sub("_", prefix) + ".csv"
        n = write_csv(records, out_dir / "tcga_subtype" / fname)
        studies = sorted(grp[STUDY_COL].dropna().astype(str).unique())
        labels = sorted(grp["_label"].unique())
        summary.append((prefix, n, studies, labels))

    print(f"\n    {'ファイル':<16}{'患者数':>7}  対象がん種 / ラベル")
    for prefix, n, studies, labels in summary:
        note = " ← PAM50" if prefix == "BRCA" else ""
        print(f"    {prefix + '.csv':<16}{n:>7}  {','.join(studies)}{note}")
        print(f"    {'':<16}{'':>7}  {', '.join(labels)}")

    print("\n使い方の例(乳がん):")
    print('    python3 preprocess_tcga.py --cancer-type BRCA ... \\')
    print(f'        --tumor-subtype-table "Immune Subtype (Thorsson 2018)='
          f'{out_dir}/Immune_Subtype.csv:2" \\')
    print(f'        --tumor-subtype-table "PAM50 (TCGA)='
          f'{out_dir}/tcga_subtype/BRCA.csv:1"')


if __name__ == "__main__":
    main()
