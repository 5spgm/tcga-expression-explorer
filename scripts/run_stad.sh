#!/usr/bin/env bash
# =====================================================================
# 胃がん(STAD)を追加する
# =====================================================================
#   bash scripts/run_stad.sh 2>&1 | tee run_stad_$(date +%Y%m%d_%H%M).log
#
# 3世代の行列は作成済みという前提。
# 新世代をSTAR-Countsから作り直す場合は先に:
#   python3 scripts/make_star_matrix.py --sample-sheet <sheet> \
#       --data-dir <GDCの展開先> --value tpm_unstranded --output <出力.tsv>
#
# 前提: run_production.sh のセクション1でサブタイプ表が生成済みであること
set -uo pipefail

# ---------------------------------------------------------------------
# 【設定】
# ---------------------------------------------------------------------
OUT_DIR="/media/nighthawk/_mnt_16TB/webtool/tcga_data"
REF_DIR="/media/nighthawk/_mnt_16TB/webtool/tcga_ref"
MAT_DIR="/media/nighthawk/_mnt_16TB/webtool/tcga_matrix"

STAD_NEW="$MAT_DIR/TCGA-STAD_TPM.xlsx"            # 新世代 TPM
STAD_MID="$MAT_DIR/20260905_STAD_FPKM_UQ.xlsx"    # 中期 FPKM-UQ(作り直した407検体版)
STAD_OLD="$MAT_DIR/TCGA-STAD-RNAseq.txt"          # 旧世代 Xena形式 (?|EntrezID)

# FFPE検体の判定に使う。GDCからダウンロードしたときのsample sheet。
# ファイル名は "gdc_sample_sheet.2026-09-05.tsv"(ドット)と
# "gdc_sample_sheet_2026-09-05.tsv"(アンダースコア)の両方がありうるので、
# 見つかった方を使う。
GDC_DIR="/home/nighthawk/X690_E/TCGA_clinical_samples/TCGA-STAD/20260905_STAD_RNASeq"
GDC_SHEET=""
for cand in "$GDC_DIR"/gdc_sample_sheet*.tsv; do
    [ -f "$cand" ] && GDC_SHEET="$cand" && break
done

# =====================================================================
# 0. 事前チェック
# =====================================================================
echo "================ 0. 事前チェック ================"
MISSING=()
check() {
    if [ -e "$2" ]; then printf "  [OK]   %-24s %s\n" "$1" "$2"
    else printf "  [無い] %-24s %s\n" "$1" "$2"; MISSING+=("$1: $2"); fi
}
check "新世代 TPM"          "$STAD_NEW"
check "中期 FPKM-UQ"        "$STAD_MID"
check "旧世代 (Xena形式)"    "$STAD_OLD"
if [ -n "$GDC_SHEET" ]; then
    printf "  [OK]   %-24s %s\n" "GDC sample sheet" "$GDC_SHEET"
else
    printf "  [無い] %-24s %s\n" "GDC sample sheet" "$GDC_DIR/gdc_sample_sheet*.tsv"
    MISSING+=("GDC sample sheet: $GDC_DIR")
fi
check "Liu 分子サブタイプ"   "$REF_DIR/subtypes/Molecular_Subtype.csv"
check "免疫サブタイプ表"     "$REF_DIR/subtypes/Immune_Subtype.csv"
check "iCluster表"          "$REF_DIR/icluster_subtypes.csv"

if [ ${#MISSING[@]} -gt 0 ]; then
    echo; echo "以下が見つかりません:"; printf '  - %s\n' "${MISSING[@]}"; exit 1
fi

# xlsxを2つ読むので、拡張子と同名の .tsv があればそちらの方が速い
for f in "$STAD_NEW" "$STAD_MID"; do
    tsv="${f%.xlsx}.tsv"
    [ -f "$tsv" ] && echo "  [ヒント] $tsv があります。指定を差し替えると読み込みが数分短縮できます。"
done
echo "  すべて揃っています。"
echo "  (xlsxは60,660行の解析に数分かかります。進まなくても異常ではありません)"
set -e

# =====================================================================
# 1. FFPE検体の除外リスト
# =====================================================================
# 2026年版のsample sheetは "Preservation Method" 列を持つ。
# vial文字(01A/01B)ではFFPEを判別できない(乳がんでは01AのFFPEが1件あった)。
# 該当0件でも空ファイルができるだけで害はない。
# 胃がんの2026-09-05版では Preservation Method が Unknown 372 / OCT 76 で、
# FFPEは0件だった。0件でも空ファイルができるだけで前処理は正常に通る。
echo; echo "================ 1. FFPE除外リスト ================"
python3 scripts/make_ffpe_exclude_list.py \
    --input "$GDC_SHEET" --output "$REF_DIR/ffpe_exclude_stad.txt" || true
FFPE_N=$(grep -c . "$REF_DIR/ffpe_exclude_stad.txt" 2>/dev/null || echo 0)
echo "  FFPE検体: $FFPE_N 件$([ "$FFPE_N" = "0" ] && echo "(除外対象なし。正常です)")"

# =====================================================================
# 2. 3世代を統合
# =====================================================================
# 胃がんは既存のサブタイプ表がそのまま5種類使える。追加の変換は不要。
#   Molecular Subtype  CIN 223 / MSI 73 / GS 50 / EBV 30 / HM-SNV 7
#   MSI Status         MSS 256 / MSI-H 75 / MSI-L 52
#   CIMP               Non-CIMP 220 / GEA CIMP-L 68 / CIMP-H 63 / CIMP EBV 31
#                      ※ CRC CIMP-L が1例だけ紛れるため min-n を 2 にして Other へ寄せる
#   Immune Subtype     C2 IFN-gamma 210 / C1 129 / C3 36
#   iCluster           C20 178 / C18:pan-GI (MSI) 115 / C1:STAD (EBV-CIMP) 32
#
# EBV陽性胃癌が CIMP EBV / C1:STAD (EBV-CIMP) として独立して現れるため、
# 大腸がんのCIMPと同じ軸で消化管のCIMPを見比べられる。
echo; echo "================ 2. STAD の前処理 ================"
python3 scripts/preprocess_tcga.py \
    --cancer-type STAD \
    --new-tpm       "$STAD_NEW" \
    --mid-fpkmuq    "$STAD_MID" \
    --old-normcount "$STAD_OLD" \
    --dedup-vials first \
    --exclude-samples "$REF_DIR/ffpe_exclude_stad.txt" \
    --tumor-subtype-table "Molecular Subtype (Liu 2018)=$REF_DIR/subtypes/Molecular_Subtype.csv:1" \
    --tumor-subtype-table "MSI Status (Liu 2018)=$REF_DIR/subtypes/MSI_Status.csv:1" \
    --tumor-subtype-table "CIMP (Liu 2018)=$REF_DIR/subtypes/CIMP.csv:2" \
    --tumor-subtype-table "Immune Subtype (Thorsson 2018)=$REF_DIR/subtypes/Immune_Subtype.csv:2" \
    --tumor-subtype-table "iCluster (Hoadley 2018)=$REF_DIR/icluster_subtypes.csv:5" \
    --entrez-cache "$REF_DIR/entrez_cache.json" \
    --out-dir "$OUT_DIR"

# =====================================================================
# 3. 圧縮 → R2用ツリー → アップロード
# =====================================================================
echo; echo "================ 3. gzip ================"
find "$OUT_DIR/STAD" -name '*.json' ! -name '*.gz' -print0 \
    | xargs -0 -P 4 -n 50 gzip -k9 -f
gzip -k9 -f "$OUT_DIR/manifest.json"

UP="/media/nighthawk/_mnt_16TB/webtool/r2_upload"
echo; echo "================ 4. R2用ツリーの更新 ================"
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
print(json.dumps(m['cancer_types']['STAD'], ensure_ascii=False, indent=2))
"
