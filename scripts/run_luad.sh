#!/usr/bin/env bash
# =====================================================================
# 肺腺がん(LUAD)を追加する
# =====================================================================
#   bash scripts/run_luad.sh 2>&1 | tee run_luad_$(date +%Y%m%d_%H%M).log
#
# 前提: run_production.sh のセクション1でサブタイプ表が生成済みであること
set -uo pipefail

# ---------------------------------------------------------------------
# 【設定】
# ---------------------------------------------------------------------
OUT_DIR="/media/nighthawk/_mnt_16TB/webtool/tcga_data"
REF_DIR="/media/nighthawk/_mnt_16TB/webtool/tcga_ref"
MAT_DIR="/media/nighthawk/_mnt_16TB/webtool/tcga_matrix"

# 新世代: GDCからダウンロードしたSTAR-Countsの置き場
GDC_DIR="/home/nighthawk/X690_E/TCGA_clinical_samples/TCGA-LUAD/20260827_RNASeq"
# sample sheet は上のディレクトリか、その近くにあるはず。実際の名前に合わせる。
GDC_SHEET="$GDC_DIR/gdc_sample_sheet.2026-08-27.tsv"

LUAD_NEW="$MAT_DIR/TCGA_LUAD_TPM.tsv"          # ステップ1で生成
LUAD_MID="$MAT_DIR/LUAD_FPKM_UQ.xls"           # 中身はTSV(拡張子が.xlsでも読める)
LUAD_OLD="$MAT_DIR/LUAD-Hiseq-v2.txt"    # Entrez ID付き(xlsx版は変換時に欠落していた)

# =====================================================================
# 0. 事前チェック
# =====================================================================
echo "================ 0. 事前チェック ================"
MISSING=()
check() {
    if [ -e "$2" ]; then printf "  [OK]   %-22s %s\n" "$1" "$2"
    else printf "  [無い] %-22s %s\n" "$1" "$2"; MISSING+=("$1: $2"); fi
}
check "GDCデータ"        "$GDC_DIR"
check "GDC sample sheet" "$GDC_SHEET"
check "中期 FPKM-UQ"     "$LUAD_MID"
check "旧世代 HiseqV2"   "$LUAD_OLD"
check "免疫サブタイプ表" "$REF_DIR/subtypes/Immune_Subtype.csv"
check "iCluster表"       "$REF_DIR/icluster_subtypes.csv"

if [ ${#MISSING[@]} -gt 0 ]; then
    echo; echo "以下が見つかりません:"; printf '  - %s\n' "${MISSING[@]}"
    echo
    echo "sample sheet の実際の名前を探すには:"
    echo "  find '$GDC_DIR' -maxdepth 2 -name '*sample_sheet*'"
    exit 1
fi
echo "  すべて揃っています。"
set -e

# =====================================================================
# 1. 新世代: STAR-Counts を1つのTPM行列にまとめる
# =====================================================================
echo; echo "================ 1. 新世代TPM行列 ================"
python3 scripts/make_star_matrix.py \
    --sample-sheet "$GDC_SHEET" \
    --data-dir     "$GDC_DIR" \
    --value        tpm_unstranded \
    --output       "$LUAD_NEW"

# =====================================================================
# 1b. FFPE検体の除外リスト
# =====================================================================
# 2026年8月版のsample sheetは新形式で "Preservation Method" 列を持つ。
# FFPE由来のRNAは品質が低く発現値が信用できないため除外する。
# vial文字(01A/01B)ではFFPEを判別できない(BRCAでは01AのFFPEが1件あった)。
echo; echo "================ 1b. FFPE除外リスト ================"
python3 scripts/make_ffpe_exclude_list.py \
    --input "$GDC_SHEET" --output "$REF_DIR/ffpe_exclude_luad.txt"
FFPE_N=$(grep -c . "$REF_DIR/ffpe_exclude_luad.txt" 2>/dev/null || echo 0)
echo "  FFPE検体: $FFPE_N 件"

# =====================================================================
# 2. 3世代を統合
# =====================================================================
# ・旧世代は Xena形式("?|100130426")。symbol不明の29行はEntrez IDから
#   mygene.info で解決される(ネット接続が必要)。
# ・中期ファイルは拡張子が .xls だが中身はTSVで、R の write.table() による
#   列数不一致(ヘッダー572 / データ573)がある。自動検出して復元される。
# ・旧世代に再発巣(02)が2検体ある。既定では除外され、サイト側の
#   チェックボックスで合流できる。
# ・LUADの既発表サブタイプ(Thorsson の TCGA Subtype 列)は 1〜6 の数値
#   のままで意味が公開されていないため、採用していない。
echo; echo "================ 2. LUAD の前処理 ================"
python3 scripts/preprocess_tcga.py \
    --cancer-type LUAD \
    --new-tpm       "$LUAD_NEW" \
    --mid-fpkmuq    "$LUAD_MID" \
    --old-normcount "$LUAD_OLD" \
    --dedup-vials first \
    --exclude-samples "$REF_DIR/ffpe_exclude_luad.txt" \
    --tumor-subtype-table "Immune Subtype (Thorsson 2018)=$REF_DIR/subtypes/Immune_Subtype.csv:2" \
    --tumor-subtype-table "iCluster (Hoadley 2018)=$REF_DIR/icluster_subtypes.csv:5" \
    --entrez-cache "$REF_DIR/entrez_cache.json" \
    --out-dir "$OUT_DIR"

# =====================================================================
# 3. 圧縮 + R2へアップロード
# =====================================================================
echo; echo "================ 3. gzip ================"
find "$OUT_DIR/LUAD" -name '*.json' ! -name '*.gz' -print0 \
    | xargs -0 -P 4 -n 50 gzip -k9 -f
gzip -k9 -f "$OUT_DIR/manifest.json"

# gzipの中身を .json という名前でR2に上げる(Content-Encodingで展開させる)
UP="/media/nighthawk/_mnt_16TB/webtool/r2_upload"
echo; echo "================ 4. R2用ツリーの更新 ================"
mkdir -p "$UP/LUAD"
( cd "$OUT_DIR" && find . -name '*.json.gz' -printf '%P\n' | while IFS= read -r f; do
      mkdir -p "$UP/$(dirname "${f%.gz}")"; ln -f "$f" "$UP/${f%.gz}"
  done )
echo "  ファイル数: $(find "$UP" -type f | wc -l)"

echo; echo "================ 5. R2へアップロード ================"
rclone copy "$UP" r2:tcga-data \
    --header-upload "Content-Encoding: gzip" \
    --header-upload "Content-Type: application/json" \
    --transfers 32 --checkers 32 --progress

echo
echo "================ 完了 ================"
python3 -c "
import json
m = json.load(open('$OUT_DIR/manifest.json'))
print('登録されているがん種:', ', '.join(sorted(m['cancer_types'])))
print()
print(json.dumps(m['cancer_types']['LUAD'], ensure_ascii=False, indent=2))
"
