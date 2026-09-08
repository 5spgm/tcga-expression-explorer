#!/usr/bin/env bash
# =====================================================================
# Zenodo寄託用パッケージを作る(全がん種まとめて1レコード)
# =====================================================================
#   bash scripts/build_zenodo.sh 2>&1 | tee build_zenodo_$(date +%Y%m%d).log
#
# がん種を増やしたら、下の MATRICES と PROVENANCE に3行ずつ足して
# VERSION を上げ、Zenodoでは「New version」として追加する。
# Concept DOI(全版を指すDOI)は変わらないので、論文の引用は壊れない。
set -uo pipefail

# ---------------------------------------------------------------------
# 【設定】
# ---------------------------------------------------------------------
VERSION="1.0"
OUT="/media/nighthawk/_mnt_16TB/webtool/zenodo_v${VERSION}"

REF="/media/nighthawk/_mnt_16TB/webtool/tcga_ref"
MAT="/media/nighthawk/_mnt_16TB/webtool/tcga_matrix"
TC="$HOME/X690_E/TCGA_clinical_samples"
REPO="./"

# --- 発現行列: "がん種:世代=ファイル名:値の種類" ------------------------
# ファイル名は $MAT からの相対パス。実際の名前に合わせて直すこと。
MATRICES=(
  "BRCA:old=20260823_BRCA-HiseqV2-Tumor-Normal.xlsx:normalized_count"
  # 中身は log2(FPKM+1)。UCSC Xena GDC hub が log2 変換して配布しているものを
  # そのまま置いている。値の種類を "FPKM" と書くと寄託ファイル名と inventory.json が
  # 中身を偽ることになるので、スケールまで名前に出す。
  # 再現時は preprocess_tcga.py に --mid-fpkm と --mid-log2p1 の両方が必要。
  # 中期の乳がんだけは UCSC Xena から取得した原本 .gz をそのまま寄託する。
  # to_tsv_gz は入力が既にgzipならbyte-for-byteでコピーするので、Excel往復の
  # 丸め(最大 3.2e-14)が入らず、gzipヘッダ内の配布側の名前と時刻も残る。
  # tcga_matrix/TCGA-BRCA.htseq_fpkm.tsv.gz は原本へのシンボリックリンク
  # (原本のファイル名にはブラウザ由来の " (1)" が付いている)。
  "BRCA:mid=TCGA-BRCA.htseq_fpkm.tsv.gz:FPKM_log2p1"
  "BRCA:new=TCGA_BRCA_TPM.xlsx:TPM"

  "COAD:old=COAD-HiseqV2-20150129.txt:normalized_count"
  "COAD:mid=COAD_FPKM_UQ.xls:FPKM_UQ"
  "COAD:new=20221216_COAD_TPM.xls:TPM"

  "PAAD:old=TCGA-PANC-RNAseq.txt:normalized_count"
  "PAAD:mid=20260824_PAAD_FPKM_UQ.xlsx:FPKM_UQ"
  "PAAD:new=20260824_TCGA_PAAD_TPM.xlsx:TPM"

  "LUAD:old=LUAD-Hiseq-v2.txt:normalized_count"
  "LUAD:mid=LUAD_FPKM_UQ.xls:FPKM_UQ"
  "LUAD:new=TCGA_LUAD_TPM.tsv:TPM"

  "STAD:old=TCGA-STAD-RNAseq.txt:normalized_count"
  "STAD:mid=20260905_STAD_FPKM_UQ.xlsx:FPKM_UQ"
  "STAD:new=TCGA-STAD_TPM.xlsx:TPM"
)

# --- 取得時の記録: "がん種:世代=ディレクトリ" ---------------------------
# sample sheet / manifest が残っているディレクトリ。
# 場所が分からないものは scripts/find_provenance.sh で探せる。
# 存在しないものは警告が出るだけで止まらないので、空欄のままでも実行できる。
# 15層のうち14層。各行列の検体集合と各シートの Sample ID を突き合わせて決めた
# (out/provenance_reconciliation.csv)。ただし照合が決定的なのは中期・新世代だけで、
# 旧世代はバーコードにvial文字が無く患者+sample typeまでしか一致を見られないため、
# 旧世代は「RNASeqV2のDCC一式が実在するディレクトリ」を根拠に置いている
# (file_manifest.txt の Platform Type が全行 RNASeqV2 / IlluminaHiSeq_RNASeqV2 であることを確認済み)。
#
# 埋まっていない層:
#   BRCA:mid  GDCから取得していない(UCSC Xena GDC hub由来)ためsample sheetが存在しない
#   PAAD:old  RNASeqV2のDCC manifestが手元に見つからない
#             ($TC/TCGA-PAAD/file_manifest.txt は Clinical(Biotab) 17行のみ)
PROVENANCE=(
  "BRCA:old=$TC/TCGA-BRCA/BRCA-HiseqV2"
  "BRCA:new=$TC/TCGA-BRCA/20250626_RNASeq"

  "COAD:old=$TC/TCGA-Colon_rectal_cancer/COAD-HiseqV2"
  "COAD:mid=$TC/TCGA-COAD/20210205_COAD_FPKM_UQ"
  "COAD:new=$TC/TCGA-COAD/20221216_COAD_STAR_counts"

  "LUAD:old=$TC/TCGA-LUAD"
  "LUAD:mid=$TC/TCGA-LUAD/20210205_LUAD_FPKM_UQ"
  "LUAD:new=$TC/TCGA-LUAD/20260827_RNASeq"

  # $TC/TCGA-PAAD/ 直下には manifest が2本あり、ディレクトリごと拾うと
  # Clinical(Biotab) の方まで旧世代の記録として寄託してしまうので、
  # RNASeqV2 の一式だけをファイル単位で指定する。
  # file_manifest (2).txt : RNASeqV2 Level3 264行 (= 44検体 x 6ファイル種)。
  #   44検体は TCGA-PANC-RNAseq.txt の44検体と完全一致する。
  #   寄託時は括弧付きの名前を避けて file_manifest.txt に付け替える
  #   (元パスは inventory.json の source_path に残る)。
  "PAAD:old=$TC/TCGA-PAAD/file_manifest (2).txt::file_manifest.txt"
  "PAAD:old=$TC/TCGA-PAAD/FILE_SAMPLE_MAP.txt"
  "PAAD:old=$TC/TCGA-PAAD/file_annotations.txt"
  "PAAD:old=$TC/TCGA-PAAD/README_DCC.txt"

  "PAAD:mid=$TC/TCGA-PAAD/20210205_PAAD_FPKM_UQ"
  "PAAD:new=$TC/TCGA-PAAD/20230221_PAAD_RNAseq"

  "STAD:old=$TC/TCGA-STAD"
  "STAD:mid=$TC/TCGA-STAD/20210205_STAD_FPKM_UQ"
  # 新世代の行列は 2024-10-14 取得分から作っている。2026-09-05 の再出力もあるが
  # 検体集合は同一で、行列の出所は 2024-10-14 の方。
  "STAD:new=$TC/TCGA-STAD/20241014_STAD_RNASeq"
)

# --- 除外リスト: "がん種=ファイル" --------------------------------------
EXCLUDES=(
  "BRCA=$REF/ffpe_exclude_brca.txt"
  "PAAD=$REF/subtypes/paad/PAAD_exclude_all.txt"
  "LUAD=$REF/ffpe_exclude_luad.txt"
  "STAD=$REF/ffpe_exclude_stad.txt"
)

# =====================================================================
# 0. 事前チェック(行列だけは無いと話にならないので厳しく見る)
# =====================================================================
echo "================ 0. 事前チェック ================"
MISSING=0
for spec in "${MATRICES[@]}"; do
    fname="${spec#*=}"; fname="${fname%:*}"
    if [ -f "$MAT/$fname" ]; then
        printf "  [OK]   %s\n" "${spec%%=*}  $fname"
    else
        printf "  [無い] %s\n" "${spec%%=*}  $MAT/$fname"
        MISSING=$((MISSING+1))
    fi
done
if [ "$MISSING" -gt 0 ]; then
    echo; echo "$MISSING 件の行列が見つかりません。MATRICES のファイル名を直してください。"
    echo "参考: ls $MAT"
    exit 1
fi
echo "  行列 ${#MATRICES[@]} 件すべて揃っています。"
set -e

# =====================================================================
# 1. パッケージを組み立てる
# =====================================================================
ARGS=()
for m in "${MATRICES[@]}";   do ARGS+=(--matrix "$m"); done
for p in "${PROVENANCE[@]}"; do ARGS+=(--provenance "$p"); done
for e in "${EXCLUDES[@]}";   do
    f="${e#*=}"; [ -f "$f" ] && ARGS+=(--exclude-list "$e")
done

echo; echo "================ 1. パッケージ作成 ================"
echo "  (xlsxはTSVに変換するため、1ファイルあたり数分かかります)"
rm -rf "$OUT"
python3 "$REPO/scripts/make_zenodo_package.py" \
    --out-dir     "$OUT" \
    --matrix-dir  "$MAT" \
    --subtype-dir "$REF" \
    --version     "$VERSION" \
    "${ARGS[@]}"

# =====================================================================
# 2. READMEを同梱する
# =====================================================================
cp "$REPO/scripts/ZENODO_README_en.md" "$OUT/README.md"
cp "$REPO/scripts/ZENODO_README_ja.md" "$OUT/README_ja.md"

# =====================================================================
# 3. アップロード用にまとめる
# =====================================================================
# Zenodoはフォルダ構造を保持しないため、本体はtar.gzにする。
# READMEは画面上で読めるよう、外に出しておく。
echo; echo "================ 2. アップロード用ファイル ================"
PARENT="$(dirname "$OUT")"
TARBALL="$PARENT/tcga-expression-explorer-data-v${VERSION}.tar.gz"
tar czf "$TARBALL" -C "$OUT" \
    matrix sample_annotation provenance inventory.json checksums.sha256 \
    $([ -d "$OUT/exclude" ] && echo exclude)
sha256sum "$TARBALL" > "$TARBALL.sha256"

echo "  $TARBALL"
echo "    $(du -h "$TARBALL" | cut -f1) / $(tar tzf "$TARBALL" | wc -l) ファイル"
echo
echo "Zenodoにアップロードするもの:"
echo "  1. $TARBALL"
echo "  2. $TARBALL.sha256"
echo "  3. $OUT/README.md"
echo "  4. $OUT/README_ja.md"
