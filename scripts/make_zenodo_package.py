#!/usr/bin/env python3
"""
make_zenodo_package.py
======================
Zenodoに寄託するデータパッケージを組み立てる。

## 設計方針

発現行列は**取得したままの形**で置く。世代ごとに遺伝子ID体系も検体集合も
違うため、1つの行列に統合するとNAだらけになり(膵臓がんの旧世代は
66%の行・76%の列がNA)、しかもID変換という「こちらの解釈」が
元データに混ざって検証できなくなる。

代わりに、解釈は別ファイルへ分ける:

    matrix/               世代ごとの発現行列(そのまま)
    sample_annotation/    検体ごとの属性・除外理由・サブタイプ割り当て
    mapping/              旧世代のID対応表(SYMBOL|Entrez -> Ensembl)
    example/              数遺伝子分の縦持ちデータ(構造の見本)
    README.md             出典・作成手順・各ファイルの説明
    checksums.sha256

## 使い方

    python3 scripts/make_zenodo_package.py \\
        --data-dir  /media/.../tcga_data \\
        --matrix-dir /media/.../tcga_matrix \\
        --subtype-dir /media/.../tcga_ref/subtypes \\
        --out-dir   /media/.../zenodo_v1.0 \\
        --matrix "BRCA:old=BRCA-HiseqV2.txt:normalized_count" \\
        --matrix "BRCA:mid=TCGA-BRCA_htseq_fpkm.tsv:FPKM" \\
        --matrix "BRCA:new=TCGA_BRCA_TPM.tsv:TPM" \\
        --example-genes TP53,KRAS,MYEOV

--matrix は「がん種:世代=ファイル名:値の種類」の形式。複数回指定する。
ファイル名は --matrix-dir からの相対パス。
"""

import argparse
import gzip
import hashlib
import json
import shutil
import sys
from datetime import datetime
from collections import defaultdict
from pathlib import Path

import pandas as pd

GENERATION_LABELS = {
    "old": "UNC IlluminaHiSeq RNASeqV2 (pre-2016)",
    "mid": "GDC HTSeq (DR15-31)",
    "new": "GDC STAR-Counts (DR32+)",
}


def sha256(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def to_tsv_gz(src: Path, dst: Path) -> None:
    """発現行列を TSV.gz として書き出す。

    寄託物は「拡張子が中身を正しく表している」ことが重要なので、
    元が xlsx / xls の場合はTSVに変換してから圧縮する。
    そのまま圧縮して .tsv.gz と名付けると、中身がxlsxなのに拡張子はTSVという
    嘘のファイルになり、受け取った人がpandasで読めずに詰まる。
    (変換には60,660行で数分かかる)
    """
    dst.parent.mkdir(parents=True, exist_ok=True)

    # 中身の判定は拡張子ではなく先頭バイトで行う。
    # このプロジェクトには「中身がTSVなのに拡張子が .xls」のファイルが実在する
    # (COAD_FPKM_UQ.xls / 20221216_COAD_TPM.xls / LUAD_FPKM_UQ.xls の3件)。
    # 拡張子で pd.read_excel に渡すと
    #   ValueError: Excel file format cannot be determined
    # で寄託パッケージのビルドが止まる。
    with src.open("rb") as fi:
        magic = fi.read(8)
    if magic[:4] == b"PK\x03\x04":
        engine = "openpyxl"          # xlsx (zip)
    elif magic[:8] == b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1":
        engine = "xlrd"              # 旧 xls (OLE2)
    else:
        engine = None                # プレーンテキスト or gzip

    # 拡張子と中身が食い違うファイルは寄託物の由来として記録に残す
    ext = src.suffix.lower()
    if (engine is None and ext in (".xlsx", ".xls")) or (
            engine is not None and ext in (".txt", ".tsv", ".csv")):
        print(f"      [注意] {src.name}: 拡張子 {ext} だが中身は "
              f"{'テキスト' if engine is None else 'Excel'}。中身に従って処理します",
              flush=True)

    if engine is not None:
        print(f"      {src.name} をTSVに変換中…(数分かかります)", flush=True)
        df = pd.read_excel(src, engine=engine)
        # R の write.table 由来で遺伝子ID列がインデックスに吸われている場合を戻す
        if df.index.name is not None or not isinstance(df.index, pd.RangeIndex):
            df = df.reset_index()
        # 名前の無い先頭列(遺伝子ID)に名前を付ける。"Unnamed: 0" のまま
        # 配ると、受け取った人が何の列か分からない。
        first = str(df.columns[0])
        if first.startswith("Unnamed") or first in ("", "index"):
            df = df.rename(columns={df.columns[0]: "gene_id"})
        with gzip.open(dst, "wt", compresslevel=9, encoding="utf-8", newline="") as fo:
            df.to_csv(fo, sep="\t", index=False)
        return

    if magic[:2] == b"\x1f\x8b":       # 既にgzip
        shutil.copyfile(src, dst)
        return
    with src.open("rb") as fi, gzip.open(dst, "wb", compresslevel=9) as fo:
        shutil.copyfileobj(fi, fo, length=1 << 20)


def parse_matrix_spec(spec: str) -> tuple[str, str, str, str]:
    """'BRCA:old=BRCA-HiseqV2.txt:normalized_count' を分解する"""
    try:
        left, right = spec.split("=", 1)
        cancer, gen = left.split(":", 1)
        fname, value_type = right.rsplit(":", 1)
    except ValueError:
        sys.exit(f"--matrix の形式が不正です: {spec}\n"
                 f"  正しい形式: がん種:世代=ファイル名:値の種類")
    if gen not in GENERATION_LABELS:
        sys.exit(f"世代は old / mid / new のいずれかです: {gen}")
    return cancer.strip(), gen.strip(), fname.strip(), value_type.strip()


# ----------------------------------------------------------------------
# 検体アノテーション表
# ----------------------------------------------------------------------
def build_sample_annotation(cancer: str, data_dir: Path, subtype_dir: Path,
                            matrices: dict) -> pd.DataFrame:
    """検体ごとに1行の表を作る。どの世代に存在するか、どう分類されたか、
    除外されたなら理由は何かを1枚で追えるようにする。
    ここが「積み上げた判断」を検証可能にする中心のファイル。

    注意: 世代によってバーコードの粒度が違う。
      旧世代(Xena)  TCGA-F2-6880-01    vial文字なし
      中期/新世代    TCGA-F2-6880-01A   vial文字あり
    そのまま突き合わせると同一検体が別行に割れるため、
    「患者ID + sample type code」を照合キーにし、実際のバーコードは
    世代ごとの列に保持する。
    """
    import re
    BC = re.compile(r"^(TCGA-\w{2}-\w{4})-(\d{2})([A-Z]?)$")
    PRIMARY = {"01", "03", "05", "09"}
    EXTRA = {"02", "04", "06", "07", "08", "40"}
    GENS = ("old", "mid", "new")

    # (患者ID, sample type code) -> {世代: バーコード}
    found: dict[tuple[str, str], dict[str, str]] = defaultdict(dict)
    for gen, path in matrices.items():
        try:
            cols = pd.read_csv(path, sep="\t", nrows=0).columns
        except Exception:
            cols = pd.read_excel(path, nrows=0).columns
        for c in cols:
            m = BC.match(str(c).strip())
            if m:
                found[(m.group(1), m.group(2))][gen] = str(c).strip()

    rows = []
    for (pid, code), per_gen in sorted(found.items()):
        group = ("Tumor" if code in PRIMARY else
                 "TumorExtra" if code in EXTRA else "Normal")
        row = {
            "patient_id": pid,
            "sample_type_code": str(code).zfill(2),
            "group": group,
            "generations": ",".join(g for g in GENS if g in per_gen),
            "n_generations": sum(1 for g in GENS if g in per_gen),
        }
        for g in GENS:
            row[f"barcode_{g}"] = per_gen.get(g, "")
        rows.append(row)
    df = pd.DataFrame(rows)

    # サブタイプ割り当てを患者単位で結合する
    for csv in sorted(subtype_dir.rglob("*.csv")):
        try:
            sub = pd.read_csv(csv)
        except Exception:
            continue
        if not {"patient_id", "subtype"} <= set(sub.columns):
            continue
        mapping = dict(zip(sub["patient_id"].astype(str), sub["subtype"]))
        col = df["patient_id"].map(mapping)
        if col.notna().sum() == 0:
            continue          # このがん種には該当しない分類
        df[f"subtype_{csv.stem}"] = col.fillna("")

    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", required=True, help="パッケージの出力先")
    ap.add_argument("--matrix-dir", required=True, help="発現行列の置き場")
    ap.add_argument("--subtype-dir", required=True, help="サブタイプCSVの置き場")
    ap.add_argument("--data-dir", help="生成済みJSONの置き場(example生成に使う)")
    ap.add_argument("--matrix", action="append", default=[], required=True,
                    metavar="がん種:世代=ファイル名:値の種類",
                    help="寄託する発現行列。複数回指定する")
    ap.add_argument("--exclude-list", action="append", default=[],
                    metavar="がん種=ファイル名", help="除外リスト(任意、複数回可)")
    ap.add_argument("--provenance", action="append", default=[],
                    metavar="がん種:世代=ディレクトリまたはファイル[::寄託名]",
                    help="取得時のGDC sample sheet / manifest が残っている"
                         "ディレクトリ。ファイルを直接指定してもよく、その場合は"
                         "そのファイルだけを寄託する(用途の違う manifest が"
                         "同居しているディレクトリで使う)。'::名前' を付けると"
                         "寄託時のファイル名を変えられる。同じ がん種:世代 に対して"
                         "複数回指定できる。"
                         "ディレクトリ指定では中の *sample_sheet*.tsv と *manifest*.txt を"
                         "取得日つきで収録する。複数回指定可")
    ap.add_argument("--example-genes", default="TP53,KRAS",
                    help="縦持ちの見本に含める遺伝子シンボル(カンマ区切り)")
    ap.add_argument("--version", default="1.0", help="パッケージのバージョン")
    args = ap.parse_args()

    out = Path(args.out_dir)
    mat_dir = Path(args.matrix_dir)
    sub_dir = Path(args.subtype_dir)
    out.mkdir(parents=True, exist_ok=True)

    # --- 1. 発現行列 ---------------------------------------------------
    print("=" * 60)
    print("1. 発現行列")
    print("=" * 60)
    specs = [parse_matrix_spec(s) for s in args.matrix]
    by_cancer: dict[str, dict[str, Path]] = defaultdict(dict)
    inventory = []
    for cancer, gen, fname, vtype in specs:
        src = mat_dir / fname
        if not src.is_file():
            sys.exit(f"見つかりません: {src}")
        dst = out / "matrix" / cancer / f"{cancer}_{gen}_{vtype}.tsv.gz"
        to_tsv_gz(src, dst)
        by_cancer[cancer][gen] = src
        # 行数・列数を数える(検証用のメタ情報として残す)
        try:
            with gzip.open(dst, "rt", encoding="utf-8") as fh:
                ncol = len(fh.readline().rstrip("\n").split("\t"))
        except Exception:  # noqa: BLE001
            ncol = 0
        inventory.append({"cancer_type": cancer, "generation": gen,
                          "value_type": vtype, "file": str(dst.relative_to(out)),
                          "source_file": fname, "n_columns": ncol,
                          "size_mb": round(dst.stat().st_size / 1e6, 1)})
        print(f"  {cancer:<6} {gen:<4} {vtype:<18} -> {dst.relative_to(out)} "
              f"({dst.stat().st_size/1e6:,.1f} MB)")

    # --- 2. 検体アノテーション ------------------------------------------
    print()
    print("=" * 60)
    print("2. 検体アノテーション")
    print("=" * 60)
    for cancer, mats in by_cancer.items():
        ann = build_sample_annotation(cancer, Path(args.data_dir or "."), sub_dir, mats)
        dst = out / "sample_annotation" / f"{cancer}_samples.tsv"
        dst.parent.mkdir(parents=True, exist_ok=True)
        # sample_type_code は "01" のような先頭ゼロ付きの文字列。
        # そのまま書くと読む側で数値の 1 に解釈されるため、
        # 文字列であることが分かる形で保存する。
        ann.to_csv(dst, sep="\t", index=False)
        # 読み込み例をREADMEに書けるよう、dtypeの指針を添える
        (dst.parent / f"{cancer}_samples.README.txt").write_text(
            "sample_type_code は先頭ゼロ付きの文字列です。読み込み例:\n"
            f"  pandas: pd.read_csv('{dst.name}', sep='\\t', dtype={{'sample_type_code': str}})\n"
            f"  R:      read.delim('{dst.name}', colClasses = c(sample_type_code = 'character'))\n")
        subcols = [c for c in ann.columns if c.startswith("subtype_")]
        print(f"  {cancer}: {len(ann)} 検体 / 分類 {len(subcols)} 種 -> "
              f"{dst.relative_to(out)}")
        print(f"      群の内訳: {dict(ann['group'].value_counts())}")

    # --- 3. 除外リスト ---------------------------------------------------
    if args.exclude_list:
        print()
        print("=" * 60)
        print("3. 除外リスト")
        print("=" * 60)
        for spec in args.exclude_list:
            cancer, fname = spec.split("=", 1)
            src = Path(fname)
            if not src.is_file():
                src = sub_dir / fname
            if not src.is_file():
                print(f"  [警告] 見つかりません: {fname}")
                continue
            # ファイル名が既にがん種で始まっている場合は接頭辞を足さない。
            # 足すと PAAD_exclude_all.txt -> PAAD_PAAD_exclude_all.txt になり、
            # READMEの再現コマンドに書いたパスと食い違う。
            stem = src.name if src.name.upper().startswith(cancer.upper() + "_") \
                   else f"{cancer}_{src.name}"
            dst = out / "exclude" / stem
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, dst)
            n = sum(1 for line in dst.read_text().splitlines() if line.strip())
            print(f"  {cancer}: {n} 件 -> {dst.relative_to(out)}")

    # --- 4. 取得時の記録 -------------------------------------------------
    # GDCは更新のたびに配布内容を差し替えるため、「いつ、どのファイルを
    # 取得したか」は後から復元できない。取得時のsample sheetとmanifestが
    # 残っていれば、それ自体がその時点のGDCの記録になる。
    # 中期(HTSeq)のデータが現在の公式から取得できないことの裏づけにもなる。
    prov_records = []
    if args.provenance:
        print()
        print("=" * 60)
        print("4. 取得時の記録 (sample sheet / manifest)")
        print("=" * 60)
        for spec in args.provenance:
            try:
                left, src_dir = spec.split("=", 1)
                cancer, gen = left.split(":", 1)
            except ValueError:
                sys.exit(f"--provenance の形式が不正です: {spec}\n"
                         f"  正しい形式: がん種:世代=ディレクトリ")
            # 右辺はディレクトリでもファイルでもよい。
            # ディレクトリを指すと下の PATTERNS に合う全ファイルを拾うが、
            # 1つのディレクトリに用途の違う manifest が同居している場合がある。
            # 実例: $TC/TCGA-PAAD/ には
            #   file_manifest.txt      Clinical(Biotab) 17行 … 旧世代RNA-seqとは無関係
            #   file_manifest (2).txt  RNASeqV2 Level3 264行 … これが旧世代の記録
            # の2本があり、ディレクトリごと拾うと前者まで
            # 「旧世代の取得記録」として寄託してしまう。
            # そのためファイル単位の指定を許し、"元パス::寄託時の名前" で
            # 名前を付け替えられるようにしている
            # (元のフルパスは inventory.json の source_path に残る)。
            rename = None
            if "::" in src_dir:
                src_dir, rename = src_dir.rsplit("::", 1)
            src = Path(src_dir)
            if src.is_file():
                found_files, use_patterns = [src], False
            elif src.is_dir():
                found_files, use_patterns = [], True
            else:
                print(f"  [警告] ファイルもディレクトリもありません: {src}")
                continue
            # GDC形式(2016年以降)と、旧TCGA DCCのData Matrix形式の両方を拾う。
            # 旧DCC形式は file_manifest.txt / FILE_SAMPLE_MAP.txt /
            # file_annotations.txt の3点セットで降ってきており、特に
            # file_annotations.txt には DCC が「解析に使うな(Item flagged DNU)」と
            # フラグを立てたaliquotの記録が入っている。現在のGDCでは同じ形では
            # 参照できないため、残す価値が高い。
            PATTERNS = (
                "*sample_sheet*.tsv", "*sample_sheet*.txt",       # GDC
                "*manifest*.txt", "*manifest*.tsv",               # GDC / 旧DCC
                "FILE_SAMPLE_MAP.txt", "file_annotations.txt",    # 旧DCC
                "README_DCC.txt",
            )
            found = found_files
            if use_patterns:
                seen = set()
                for pat in PATTERNS:
                    for f in sorted(src.glob(pat)):
                        if f.name not in seen:
                            seen.add(f.name)
                            found.append(f)
            if not found:
                print(f"  [警告] {src} に sample sheet / manifest が見つかりません")
                continue
            for f in found:
                name = rename if (rename and len(found) == 1) else f.name
                dst = out / "provenance" / cancer.strip() / gen.strip() / name
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(f, dst)
                mtime = datetime.fromtimestamp(f.stat().st_mtime).strftime("%Y-%m-%d")
                n_rows = max(0, sum(1 for _ in f.open(errors="replace")) - 1)
                prov_records.append({
                    "cancer_type": cancer.strip(), "generation": gen.strip(),
                    "file": str(dst.relative_to(out)), "file_mtime": mtime,
                    "n_rows": n_rows, "source_path": str(f),
                })
                # 表示は寄託後の名前。付け替えた場合は元の名前も出す
                shown = name if name == f.name else f"{name} (元: {f.name})"
                print(f"  {cancer.strip():<6} {gen.strip():<4} {shown:<44} "
                      f"更新日 {mtime} / {n_rows} 行")

    # --- 5. 目録とチェックサム -------------------------------------------
    (out / "inventory.json").write_text(
        json.dumps({"version": args.version, "matrices": inventory,
                    "provenance": prov_records},
                   ensure_ascii=False, indent=2))

    print()
    print("=" * 60)
    print("6. チェックサム")
    print("=" * 60)
    lines = []
    for p in sorted(out.rglob("*")):
        if p.is_file() and p.name != "checksums.sha256":
            lines.append(f"{sha256(p)}  {p.relative_to(out)}")
    (out / "checksums.sha256").write_text("\n".join(lines) + "\n")
    print(f"  {len(lines)} ファイル -> checksums.sha256")

    total = sum(p.stat().st_size for p in out.rglob("*") if p.is_file())
    print()
    print(f"パッケージ: {out}  ({total/1e9:,.2f} GB / {len(lines)} ファイル)")
    print("README.md を書いてから、Zenodoにアップロードしてください。")


if __name__ == "__main__":
    main()
