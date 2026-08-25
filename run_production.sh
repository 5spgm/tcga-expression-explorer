#!/usr/bin/env bash
# =====================================================================
# TCGA Expression Explorer — 本番データ生成スクリプト
# =====================================================================
#   bash run_production.sh 2>&1 | tee run_$(date +%Y%m%d_%H%M).log
#
# ログは必ず残すこと。除外件数や重複vialの内訳は画面出力にしか出ないため、
# あとから「どういう条件で作ったデータか」を確認する唯一の手段になる。
#
# 実行すると、まず【0. 事前チェック】で必要なファイルの有無を全部まとめて
# 確認する。足りないものがあれば、その一覧を出して何もせずに終了する。
set -uo pipefail

# ---------------------------------------------------------------------
# 【設定1】ディレクトリ
# ---------------------------------------------------------------------
OUT_DIR="/media/nighthawk/_mnt_16TB/webtool/tcga_data"   # nginxが配信するデータ置き場
REF_DIR="/media/nighthawk/_mnt_16TB/webtool/tcga_ref"    # 論文・臨床情報
MAT_DIR="/media/nighthawk/_mnt_16TB/webtool/tcga_matrix" # 発現行列

# ---------------------------------------------------------------------
# 【設定2】入力ファイル名
#   実際のファイル名に合わせて書き換える。
#   Liu 2018 と Thorsson 2018 の補足資料は**どちらも 'mmc2.xlsx'** という
#   名前で配布されているので、取り違えに注意(区別できる名前にしておくと安全)。
# ---------------------------------------------------------------------
ICLUSTER_XLSX="$REF_DIR/mmc6.xlsx"                                    # Hoadley 2018
LIU_XLSX="$REF_DIR/mmc2_liu.xlsx"                                     # Liu 2018
THORSSON_XLSX="$REF_DIR/1-s2_0-S1074761318301213-mmc2.xlsx"           # Thorsson 2018 (Immunity)
BRCA_CLINICAL="$REF_DIR/nationwidechildrens_org_clinical_patient_brca.txt"
BRCA_SHEET="$REF_DIR/gdc_sample_sheet_2025-06-26.tsv"

BRCA_NEW="$MAT_DIR/TCGA_BRCA_TPM.xlsx"
BRCA_MID="$MAT_DIR/TCGA-BRCA_htseq_fpkm.xlsx"      # BRCAは FPKM (COADは FPKM-UQ)
BRCA_OLD="$MAT_DIR/BRCA-HiseqV2-Tumor-Normal.xlsx"

COAD_NEW="$MAT_DIR/20221216_COAD_TPM.xls"
COAD_MID="$MAT_DIR/COAD_FPKM_UQ.xls"
COAD_OLD="$MAT_DIR/COAD-HiseqV2-20150129.txt"

# 生成したくないがん種は "no" にする
RUN_BRCA="yes"
RUN_COAD="yes"

# =====================================================================
# 0. 事前チェック — 足りないファイルを全部まとめて報告する
# =====================================================================
echo "================ 0. 事前チェック ================"
MISSING=()
check() {   # check <説明> <パス>
    if [ -f "$2" ]; then
        printf "  [OK]   %-28s %s\n" "$1" "$2"
    else
        printf "  [無い] %-28s %s\n" "$1" "$2"
        MISSING+=("$1: $2")
    fi
}

check "iCluster (Hoadley)"  "$ICLUSTER_XLSX"
check "Liu 2018"            "$LIU_XLSX"
check "Thorsson 2018"       "$THORSSON_XLSX"
if [ "$RUN_BRCA" = "yes" ]; then
    check "BRCA 臨床(BCR)"   "$BRCA_CLINICAL"
    check "BRCA sample sheet" "$BRCA_SHEET"
    check "BRCA 新世代 TPM"   "$BRCA_NEW"
    check "BRCA 中期 FPKM"    "$BRCA_MID"
    check "BRCA 旧世代"       "$BRCA_OLD"
fi
if [ "$RUN_COAD" = "yes" ]; then
    check "COAD 新世代 TPM"   "$COAD_NEW"
    check "COAD 中期 FPKM-UQ" "$COAD_MID"
    check "COAD 旧世代"       "$COAD_OLD"
fi

if [ ${#MISSING[@]} -gt 0 ]; then
    echo
    echo "以下のファイルが見つかりません。パスを直すか、ファイルを配置してください。"
    printf '  - %s\n' "${MISSING[@]}"
    echo
    echo "参考: そのディレクトリにある .xlsx / .xls / .tsv / .txt"
    ls -1 "$REF_DIR" "$MAT_DIR" 2>/dev/null | grep -Ei '\.(xlsx|xls|tsv|txt)$' | sed 's/^/    /'
    exit 1
fi

echo "  すべて揃っています。"
set -e   # ここから先はエラーで即停止する
mkdir -p "$OUT_DIR" "$REF_DIR/subtypes"

# =====================================================================
# 1. サブタイプ表を作る(全がん種で使い回すので1回だけ)
# =====================================================================
echo; echo "================ 1. サブタイプ表 ================"
python3 make_icluster_subtype_table.py \
    --input "$ICLUSTER_XLSX" --output "$REF_DIR/icluster_subtypes.csv"

python3 make_liu_subtype_tables.py \
    --input "$LIU_XLSX" --out-dir "$REF_DIR/subtypes"

# 免疫サブタイプ C1-C6(全がん種共通)と、各がん種の既発表サブタイプ
# (乳がんは PAM50)を一度に生成する
python3 make_immune_subtype_tables.py \
    --input "$THORSSON_XLSX" --out-dir "$REF_DIR/subtypes"

if [ "$RUN_BRCA" = "yes" ]; then
    python3 make_tnbc_subtype_table.py \
        --input      "$BRCA_CLINICAL" \
        --output     "$REF_DIR/subtypes/BRCA_TNBC.csv" \
        --detail-csv "$REF_DIR/subtypes/BRCA_TNBC_detail.csv"
fi

# =====================================================================
# 2. FFPE検体の除外リスト
# =====================================================================
if [ "$RUN_BRCA" = "yes" ]; then
    echo; echo "================ 2. FFPE除外リスト ================"
    # GDC sample sheet の "Preservation Method" 列から作る。
    # BRCAでは13件。うち1件は 01A なので vial文字では絶対に拾えない。
    python3 make_ffpe_exclude_list.py \
        --input "$BRCA_SHEET" --output "$REF_DIR/ffpe_exclude_brca.txt"
fi

# =====================================================================
# 3. BRCA
# =====================================================================
if [ "$RUN_BRCA" = "yes" ]; then
    echo; echo "================ 3. BRCA ================"
    python3 preprocess_tcga.py \
        --cancer-type BRCA \
        --new-tpm       "$BRCA_NEW" \
        --mid-fpkm      "$BRCA_MID" \
        --old-normcount "$BRCA_OLD" \
        --exclude-samples "$REF_DIR/ffpe_exclude_brca.txt" \
        --dedup-vials first \
        --tumor-subtype-table "PAM50 (TCGA)=$REF_DIR/subtypes/tcga_subtype/BRCA.csv:1" \
        --tumor-subtype-table "TNBC (IHC)=$REF_DIR/subtypes/BRCA_TNBC.csv:1" \
        --tumor-subtype-table "Immune Subtype (Thorsson 2018)=$REF_DIR/subtypes/Immune_Subtype.csv:2" \
        --tumor-subtype-table "iCluster (Hoadley 2018)=$REF_DIR/icluster_subtypes.csv" \
        --entrez-cache "$REF_DIR/entrez_cache.json" \
        --out-dir "$OUT_DIR"
fi

# =====================================================================
# 4. COAD(転移巣の別枠化と免疫サブタイプ追加のため作り直す)
# =====================================================================
# ※ Thorsson の GI.csv (CIN/MSI/GS/EBV) は Liu の Molecular Subtype と
#   ほぼ同内容なので、両方入れるとプルダウンが重複する。Liu側を採用する。
if [ "$RUN_COAD" = "yes" ]; then
    echo; echo "================ 4. COAD ================"
    python3 preprocess_tcga.py \
        --cancer-type COAD \
        --new-tpm       "$COAD_NEW" \
        --mid-fpkmuq    "$COAD_MID" \
        --old-normcount "$COAD_OLD" \
        --dedup-vials first \
        --tumor-subtype-table "iCluster (Hoadley 2018)=$REF_DIR/icluster_subtypes.csv" \
        --tumor-subtype-table "Molecular Subtype (Liu 2018)=$REF_DIR/subtypes/Molecular_Subtype.csv:1" \
        --tumor-subtype-table "MSI Status (Liu 2018)=$REF_DIR/subtypes/MSI_Status.csv:1" \
        --tumor-subtype-table "CIMP (Liu 2018)=$REF_DIR/subtypes/CIMP.csv:1" \
        --tumor-subtype-table "Colorectal CMS (Liu 2018)=$REF_DIR/subtypes/Colorectal_CMS.csv:2" \
        --tumor-subtype-table "Immune Subtype (Thorsson 2018)=$REF_DIR/subtypes/Immune_Subtype.csv:2" \
        --entrez-cache "$REF_DIR/entrez_cache.json" \
        --out-dir "$OUT_DIR"
fi

echo
echo "================ 生成完了 ================"
cat "$OUT_DIR/manifest.json"
