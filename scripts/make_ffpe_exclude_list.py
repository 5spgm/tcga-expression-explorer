#!/usr/bin/env python3
"""
make_ffpe_exclude_list.py
=========================
FFPE由来の検体バーコード一覧を作り、preprocess_tcga.py の --exclude-samples に
渡せる形(1行1バーコード)で出力する。

FFPEかどうかはTCGAバーコードからは判別できない(vial文字 A/B/C は同じ検体の
何本目のバイアルかを示すだけで、固定法とは無関係)。判定できるのは検体レベルの
メタデータだけなので、以下のいずれかのファイルが必要になる。

**このスクリプトが受け付ける入力**(列名を自動検出する)

1. BCR biospecimen sample ファイル(推奨)
   `nationwidechildrens_org_biospecimen_sample_brca.txt`
   - バーコード列: bcr_sample_barcode
   - FFPE列      : is_ffpe (YES / NO)
   - clinical_patient ファイルには入っていないので注意。同じBCRアーカイブの
     biospecimen 側を探すこと。

2. GDC sample sheet / GDC APIの出力をCSV/TSVにしたもの
   - バーコード列: Sample ID / sample_submitter_id / submitter_id など
   - FFPE列      : is_ffpe / Is FFPE など

使い方:

    python3 make_ffpe_exclude_list.py \
        --input nationwidechildrens_org_biospecimen_sample_brca.txt \
        --output ffpe_exclude_brca.txt

    python3 preprocess_tcga.py --cancer-type BRCA ... \
        --exclude-samples ffpe_exclude_brca.txt

出力されるバーコードは `TCGA-XX-XXXX-01A` の形(sample+vialまで)に切り詰める。
preprocess_tcga.py 側の列名がこの粒度のため。
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

BARCODE_CANDIDATES = [
    "bcr_sample_barcode", "sample_id", "sample id", "sample_submitter_id",
    "submitter_id", "barcode", "bcr_aliquot_barcode",
]
FFPE_CANDIDATES = ["is_ffpe", "is ffpe", "ffpe", "is_derived_from_ffpe",
                   "preservation method", "preservation_method"]
# is_ffpe 列は YES/NO、GDC sample sheet の "Preservation Method" 列は
# FFPE / OCT / Unknown という値を取る。どちらでも拾えるようにする。
TRUTHY = {"yes", "true", "1", "y", "ffpe"}

# TCGA-XX-XXXX-01A まで(sample type 2桁 + vial文字)を取り出す
BARCODE_TRIM = re.compile(r"^(TCGA-\w{2}-\w{4}-\d{2}[A-Z]?)")


def find_column(columns, candidates, what):
    lower = {str(c).strip().lower(): c for c in columns}
    for cand in candidates:
        if cand in lower:
            return lower[cand]
    sys.exit(f"{what} に相当する列が見つかりません。\n"
             f"  探した列名: {candidates}\n"
             f"  実際の列名: {list(columns)[:25]}\n"
             f"--barcode-column / --ffpe-column で明示的に指定してください。")


def read_table(path: str) -> pd.DataFrame:
    """BCRのファイルは列名の下に説明行(別名行・CDE_ID行)が入ることがあるが、
    その行数はファイルによって違う(clinicalは2行、biospecimen sampleは1行)。
    決め打ちで読み飛ばすと先頭の実データを落とすので、CDE_IDを含む行と
    完全な空行だけを個別に検出して飛ばす。"""
    sep = "\t" if Path(path).suffix.lower() in (".txt", ".tsv") else ","
    probe = pd.read_csv(path, sep=sep, nrows=5, dtype=str, header=0, low_memory=False)
    skip = []
    for i in range(len(probe)):
        cells = [str(v) for v in probe.iloc[i].tolist()]
        joined = " ".join(cells)
        is_cde = "CDE_ID" in joined
        is_blank = all(c.strip() in ("", "nan") for c in cells)
        if is_cde or is_blank:
            skip.append(i + 1)   # ヘッダー行の次を1として数える
        else:
            break                 # 実データに到達したら終了
    if skip:
        print(f"  (説明行 {len(skip)} 行を読み飛ばしました)")
    return pd.read_csv(path, sep=sep, skiprows=skip, dtype=str, low_memory=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True,
                    help="BCR biospecimen sample ファイル、またはGDC sample sheet")
    ap.add_argument("--output", default="ffpe_exclude.txt", help="出力先")
    ap.add_argument("--barcode-column", help="バーコード列を明示指定する場合")
    ap.add_argument("--ffpe-column", help="FFPE列を明示指定する場合")
    args = ap.parse_args()

    df = read_table(args.input)
    bc_col = args.barcode_column or find_column(df.columns, BARCODE_CANDIDATES, "バーコード列")
    ffpe_col = args.ffpe_column or find_column(df.columns, FFPE_CANDIDATES, "FFPE列")
    print(f"バーコード列: '{bc_col}' / FFPE列: '{ffpe_col}'")
    print("FFPE列の値の内訳:")
    for val, n in df[ffpe_col].fillna("(空欄)").value_counts().items():
        print(f"    {val}: {n}")

    is_ffpe = df[ffpe_col].fillna("").astype(str).str.strip().str.lower().isin(TRUTHY)
    barcodes = []
    for raw in df.loc[is_ffpe, bc_col].dropna().astype(str):
        m = BARCODE_TRIM.match(raw.strip())
        if m:
            barcodes.append(m.group(1))
    barcodes = sorted(set(barcodes))

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(barcodes) + ("\n" if barcodes else ""))

    print(f"\nFFPE検体: {len(barcodes)} 件 -> {out}")
    if not barcodes:
        print("[警告] 該当が0件でした。FFPE列の値が YES/NO 以外の表記でないか、"
              "上の内訳を確認してください。")
    else:
        print(f"  例: {', '.join(barcodes[:5])}")
        print(f"\n次のように渡してください:\n"
              f"    python3 preprocess_tcga.py ... --exclude-samples {out}")


if __name__ == "__main__":
    main()
