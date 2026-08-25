#!/usr/bin/env python3
"""
check_rnaseq_list.py
====================
GDC sample sheet を読み、`20250913_RNAseq_list.R` の選別ロジックで
**どの検体・どの患者が落ちるか**を洗い出し、修正版の検体リストを出力する。

Rスクリプト側には次の3つの落とし穴がある。このスクリプトはそれぞれが実際に
何件効いているかを数え、落ちた検体を一覧表示する。

(1) `subset(sheet, sheet[,8]=="Tumor")` / `=="Normal"`
    Sample Type列を**完全一致**で見ている。GDC sample sheetの標準的な値は
    "Primary Tumor" / "Solid Tissue Normal" / "Metastatic" などなので、
    値を書き換えていない限り一致せず、一致しない行は警告なしに消える。
    特に "Metastatic" は Tumor にも Normal にも一致しないので必ず落ちる。

(2) `sheet5 <- sheet5[grep("-01A", sheet5[,7]), ]`
    コメントは「01Aと01Bが両方ある場合、01Aを落とす」だが、コードは逆に
    **01Aだけを残して他を全部落としている**。さらに、複数ファイルを持つ
    患者(sheet2系統)で 01A が1つも無い場合、その患者は**まるごと消える**。
    ファイルが1つだけの患者(sheet1)はこのgrepを通らないので、01Bしか
    無くても残る — この非対称性が発見を難しくしている。

(3) `sheet3 <- sheet3[!duplicated(sheet3[,6]), ]`
    6列目は Case ID(患者)なので、**患者単位**で重複除去している。
    同じ患者が別々の検体(01Aと06Aなど)をそれぞれ複数ファイル持つ場合、
    患者につき1行しか残らず、他の検体が消える。
    (意図はおそらく Sample ID 単位=7列目の重複除去)

使い方:

    python3 check_rnaseq_list.py \
        --sample-sheet gdc_sample_sheet.2025-09-13.tsv \
        --output corrected_file_list.tsv

出力される `corrected_file_list.tsv` は、検体(Sample ID)ごとに1ファイルを
選んだ正しい一覧。R側の `sheet` の代わりに読み込めば作り直せる。
"""

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

BARCODE_RE = re.compile(r"^(TCGA-\w{2}-\w{4}-(\d{2})([A-Z]?))")

COL_CANDIDATES = {
    "file_name": ["file name", "file_name", "filename"],
    "case_id": ["case id", "case_id", "case.id"],
    "sample_id": ["sample id", "sample_id", "sample.id"],
    # 新しいGDC sample sheetは "Tissue Type"(Tumor/Normal)、
    # 古い形式は "Sample Type"(Primary Tumor など)。両方に対応する。
    "sample_type": ["tissue type", "sample type", "sample_type", "sample.type",
                    "tissue.type"],
}


def find_col(columns, key):
    lower = {str(c).strip().lower(): c for c in columns}
    for cand in COL_CANDIDATES[key]:
        if cand in lower:
            return lower[cand]
    sys.exit(f"'{key}' に相当する列が見つかりません。実際の列名: {list(columns)}")


def parse_barcode(sample_id):
    m = BARCODE_RE.match(str(sample_id).strip())
    return (m.group(1), m.group(2), m.group(3)) if m else (None, None, None)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-sheet", required=True, help="GDC sample sheet (TSV)")
    ap.add_argument("--output", default="corrected_file_list.tsv",
                    help="修正版の検体リスト出力先")
    ap.add_argument("--prefer-vial", default="first", choices=["first", "last"],
                    help="同一検体に複数vialがある場合にどちらを採るか(既定: first=A)")
    args = ap.parse_args()

    sep = "\t" if Path(args.sample_sheet).suffix.lower() in (".tsv", ".txt") else ","
    df = pd.read_csv(args.sample_sheet, sep=sep, dtype=str).fillna("")
    c_file = find_col(df.columns, "file_name")
    c_case = find_col(df.columns, "case_id")
    c_samp = find_col(df.columns, "sample_id")
    c_type = find_col(df.columns, "sample_type")
    print(f"列: file='{c_file}' case='{c_case}' sample='{c_samp}' type='{c_type}'")
    print(f"総行数(ファイル数): {len(df)}\n")

    # ---- (1) Sample Type の完全一致でどれだけ落ちるか -----------------
    print("=" * 62)
    print("(1) Sample Type 列の値と、Rの完全一致 =='Tumor'/=='Normal' の結果")
    print("=" * 62)
    matched_rows = 0
    for val, n in df[c_type].value_counts().items():
        if val in ("Tumor", "Normal"):
            verdict = "一致 -> 残る"
            matched_rows += n
        else:
            verdict = "不一致 -> **落ちる**"
        print(f"    {val!r:<32} {n:>5} 行   {verdict}")
    dropped_by_type = len(df) - matched_rows
    if dropped_by_type:
        print(f"\n  [!] Sample Type が一致せず落ちる行: {dropped_by_type} 行")
        print("      GDC標準の値のままなら全滅する。Rを実行する前に、この列を")
        print("      'Tumor'/'Normal' に書き換えていたかどうか確認すること。")
    else:
        print("\n  この列は 'Tumor'/'Normal' に整形済み。ここでの取りこぼしは無い。")

    # 以降の解析は「Tumor扱いの行」に対して行う。列が整形済みでない場合は
    # 腫瘍系のsample type codeから推定する(検証を続けられるようにするため)。
    if matched_rows > 0:
        tumor = df[df[c_type] == "Tumor"].copy()
        normal = df[df[c_type] == "Normal"].copy()
        note = ""
    else:
        codes = df[c_samp].map(lambda s: parse_barcode(s)[1])
        tumor = df[codes.astype(str) < "10"].copy()
        normal = df[codes.astype(str) >= "10"].copy()
        note = "(Sample Type列が非標準のため、バーコードのcodeから推定した)"
    print(f"\n  Tumor行 {len(tumor)} / Normal行 {len(normal)} {note}")

    # ---- Rのロジックを再現 --------------------------------------------
    counts_case = tumor[c_case].value_counts()
    sheet1 = tumor[tumor[c_case].map(counts_case) == 1]          # 単一ファイルの患者
    sheet2 = tumor[tumor[c_case].map(counts_case) > 1]           # 複数ファイルの患者

    counts_samp = sheet2[c_samp].value_counts()
    sheet3 = sheet2[sheet2[c_samp].map(counts_samp) > 1]
    sheet3 = sheet3.drop_duplicates(subset=[c_case], keep="first")  # (3) 患者単位で重複除去
    sheet4 = sheet2[sheet2[c_samp].map(counts_samp) == 1]
    sheet5 = pd.concat([sheet3, sheet4])
    sheet5_after = sheet5[sheet5[c_samp].str.contains("-01A", na=False)]  # (2) grep

    r_selected = pd.concat([normal, sheet1, sheet5_after])
    r_samples = set(r_selected[c_samp])

    # ---- (2)(3) の影響を数える -----------------------------------------
    print()
    print("=" * 62)
    print("(2) grep('-01A') による取りこぼし")
    print("=" * 62)
    dropped_grep = sheet5[~sheet5[c_samp].str.contains("-01A", na=False)]
    print(f"    grepで落ちた検体: {len(dropped_grep)} 件")
    lost_cases = sorted(set(dropped_grep[c_case]) - set(r_selected[c_case]))
    if len(dropped_grep):
        for _, r in dropped_grep.head(20).iterrows():
            gone = " <<< この患者は他に検体が無く、まるごと消える" \
                if r[c_case] in lost_cases else ""
            print(f"      {r[c_samp]:<20} ({r[c_type]}){gone}")
        if len(dropped_grep) > 20:
            print(f"      ... 他 {len(dropped_grep) - 20} 件")
    print(f"\n  [!] 検体ごと消えた患者: {len(lost_cases)} 名")
    for c in lost_cases:
        print(f"      {c}")

    print()
    print("=" * 62)
    print("(3) 患者単位の重複除去(!duplicated(Case ID))による取りこぼし")
    print("=" * 62)
    dup_rows = sheet2[sheet2[c_samp].map(counts_samp) > 1]
    collapsed = len(dup_rows.drop_duplicates(subset=[c_samp])) - len(sheet3)
    print(f"    同一検体に複数ファイルがある行: {len(dup_rows)} 行")
    print(f"    本来は検体単位で {len(dup_rows.drop_duplicates(subset=[c_samp]))} 検体残るはずが、")
    print(f"    患者単位の除去により {len(sheet3)} 検体しか残らない -> 差分 {collapsed} 検体")

    # ---- 正しい選び方 ---------------------------------------------------
    correct_rows = []
    for sid, grp in df.groupby(c_samp, sort=False):
        correct_rows.append(grp.iloc[0])       # 1検体につきファイル1つ
    correct = pd.DataFrame(correct_rows)

    # 同一患者×同一sample typeで複数vialがあれば1つに絞る
    keep = []
    by_key = defaultdict(list)
    for _, r in correct.iterrows():
        full, code, vial = parse_barcode(r[c_samp])
        by_key[(str(r[c_case]), code)].append((vial or "", r))
    for key, items in by_key.items():
        items.sort(key=lambda t: t[0])
        keep.append(items[0][1] if args.prefer_vial == "first" else items[-1][1])
    correct = pd.DataFrame(keep)

    missing = correct[~correct[c_samp].isin(r_samples)]
    print()
    print("=" * 62)
    print("まとめ")
    print("=" * 62)
    print(f"    Rスクリプトが選んだ検体      : {len(r_samples)}")
    print(f"    本来選ばれるべき検体          : {len(correct)}")
    print(f"    取りこぼした検体              : {len(missing)}")
    lost_patients = sorted(set(missing[c_case]) - set(r_selected[c_case]))
    print(f"    完全に消えた患者              : {len(lost_patients)}")
    if len(missing):
        print("\n    取りこぼした検体の内訳(sample type別):")
        for val, n in missing[c_type].value_counts().items():
            print(f"      {val}: {n}")

    out = Path(args.output)
    correct.to_csv(out, sep="\t", index=False)
    print(f"\n修正版の検体リスト -> {out} ({len(correct)} 検体)")
    print("R側の `sheet` をこのファイルに差し替えれば、取りこぼしなく作り直せる。")


if __name__ == "__main__":
    main()
