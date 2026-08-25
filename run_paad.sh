#!/usr/bin/env bash
# =====================================================================
# 膵臓がん(PAAD)を追加する
# =====================================================================
#   bash run_paad.sh 2>&1 | tee run_paad_$(date +%Y%m%d_%H%M).log
#
# 新世代のTPM行列は作成済みという前提(make_star_matrix.py の出力)。
# 前提: run_production.sh のセクション1でサブタイプ表が生成済みであること
set -uo pipefail

# ---------------------------------------------------------------------
# 【設定】
# ---------------------------------------------------------------------
OUT_DIR="/media/nighthawk/_mnt_16TB/webtool/tcga_data"
REF_DIR="/media/nighthawk/_mnt_16TB/webtool/tcga_ref"
MAT_DIR="/media/nighthawk/_mnt_16TB/webtool/tcga_matrix"

PAAD_NEW="$MAT_DIR/TCGA_PAAD_TPM.tsv"          # 新世代(make_star_matrix.pyの出力)
PAAD_MID="$MAT_DIR/PAAD_FPKM_UQ.xlsx"          # 中期
PAAD_OLD="$MAT_DIR/TCGA-PANC-RNAseq.txt"       # 旧世代(2015年取得, Xena形式)

# TCGA-PAAD marker paper (Raphael et al., Cancer Cell 2017) の Table S1
PAAD_SUPPL="$REF_DIR/mmc2_paad.xlsx"
SUB="$REF_DIR/subtypes/paad"

# =====================================================================
# 0. 事前チェック
# =====================================================================
echo "================ 0. 事前チェック ================"
MISSING=()
check() {
    if [ -e "$2" ]; then printf "  [OK]   %-24s %s\n" "$1" "$2"
    else printf "  [無い] %-24s %s\n" "$1" "$2"; MISSING+=("$1: $2"); fi
}
check "新世代 TPM"             "$PAAD_NEW"
check "中期 FPKM-UQ"           "$PAAD_MID"
check "旧世代 2015 (Xena形式)" "$PAAD_OLD"
check "PAAD補足表 (Table S1)"  "$PAAD_SUPPL"
check "免疫サブタイプ表"       "$REF_DIR/subtypes/Immune_Subtype.csv"
if [ ${#MISSING[@]} -gt 0 ]; then
    echo; echo "以下が見つかりません:"; printf '  - %s\n' "${MISSING[@]}"; exit 1
fi
echo "  すべて揃っています。"
set -e

# =====================================================================
# 1. 膵臓がん固有のサブタイプ表と除外リストを作る
# =====================================================================
# Table S1 のシート構成:
#   FreezeSamples              解析対象150検体(分類の割り当てはここ)
#   ExcludedSamples            膵管腺癌でないと判定された19検体
#   PseudoNormals (low purity) 腫瘍細胞含有率が極端に低い9検体
#   "Real" Normals (ID -11A)   隣接正常組織4検体
# 150 + 19 + 9 = 178 で、新世代の原発巣178検体とちょうど一致する。
echo; echo "================ 1. サブタイプ表 ================"
python3 make_paad_subtype_tables.py --input "$PAAD_SUPPL" --out-dir "$SUB"

# =====================================================================
# 2. 3世代を統合して遺伝子ごとのJSONを書き出す
# =====================================================================
# --exclude-samples に PAAD_exclude_all.txt (非PDAC 19 + pseudonormal 9 = 28)
# を渡す。これにより「分類なし」表示の腫瘍数と、各分類の合計が
# どちらも150で一致し、プルダウンを切り替えてもnが動かない。
#
# 分類は4つ:
#   Moffitt      Basal-like / Classical。予後との関連が最も確立。
#                純度との交絡が少ない(high 31/45, low 34/40)。
#   Bailey       4群。ただし ADEX と Immunogenic はほぼ低純度検体で構成され
#                (ADEX 9/29, Immunogenic 2/26)、腺房細胞や免疫細胞の混入を
#                見ている可能性が指摘されている。下のPurityと併せて解釈すること。
#   Tumor purity high 76 / low 74。上記の交絡を確認するためのもの。
#   Immune       Thorsson 2018。他がん種と横断比較できる。
#
# iCluster は膵臓がんだと176検体中157件が1クラスタに集中するため入れない。
echo; echo "================ 2. PAAD の前処理 ================"
python3 preprocess_tcga.py \
    --cancer-type PAAD \
    --new-tpm       "$PAAD_NEW" \
    --mid-fpkmuq    "$PAAD_MID" \
    --old-normcount "$PAAD_OLD" \
    --dedup-vials first \
    --exclude-samples "$SUB/PAAD_exclude_all.txt" \
    --tumor-subtype-table "Moffitt (Cancer Cell 2017)=$SUB/PAAD_Moffitt.csv:1" \
    --tumor-subtype-table "Bailey (Cancer Cell 2017)=$SUB/PAAD_Bailey.csv:2" \
    --tumor-subtype-table "Tumor purity (Cancer Cell 2017)=$SUB/PAAD_Purity_Class.csv:1" \
    --tumor-subtype-table "Immune Subtype (Thorsson 2018)=$REF_DIR/subtypes/Immune_Subtype.csv:2" \
    --entrez-cache "$REF_DIR/entrez_cache.json" \
    --out-dir "$OUT_DIR"

# =====================================================================
# 3. 圧縮(gzip_static 用)
# =====================================================================
echo; echo "================ 3. gzip ================"
find "$OUT_DIR/PAAD" -name '*.json' ! -name '*.gz' -print0 \
    | xargs -0 -P 4 -n 50 gzip -k9 -f
gzip -k9 -f "$OUT_DIR/manifest.json"
echo "  完了: $(find "$OUT_DIR/PAAD" -name '*.json.gz' | wc -l) ファイル"

echo
echo "================ 生成完了 ================"
python3 -c "
import json
m = json.load(open('$OUT_DIR/manifest.json'))
print('登録されているがん種:', ', '.join(sorted(m['cancer_types'])))
print()
print(json.dumps(m['cancer_types']['PAAD'], ensure_ascii=False, indent=2))
"
