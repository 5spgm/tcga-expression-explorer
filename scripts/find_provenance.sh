#!/usr/bin/env bash
# =====================================================================
# 取得時の記録(sample sheet / manifest)を探して、取得日を推定する
# =====================================================================
#   bash scripts/find_provenance.sh ~/X690_E/TCGA_clinical_samples
#
# GDC形式と旧TCGA DCC形式の両方を探す。
# 旧DCC形式は file_annotations.txt に記録された日付の最新値が
# 「ダウンロードはそれ以降」という下限になる。
set -uo pipefail
ROOT="${1:-$HOME}"

echo "=== GDC形式 (2016年以降) ==="
find "$ROOT" -name 'gdc_sample_sheet*' -o -name 'gdc_manifest*' 2>/dev/null | sort | while read -r f; do
    printf "%s  %6s行  %s\n" \
        "$(date -r "$f" +%Y-%m-%d)" "$(wc -l < "$f")" "$f"
done

echo
echo "=== 旧TCGA DCC形式 (Data Matrix) ==="
find "$ROOT" -name 'file_manifest.txt' 2>/dev/null | sort | while read -r f; do
    d=$(dirname "$f")
    ann="$d/file_annotations.txt"
    # アノテーションの最新日付 = ダウンロード日の下限
    latest="-"
    if [ -f "$ann" ]; then
        latest=$(grep -oE '[0-9]{2}/[0-9]{2}/[0-9]{4}' "$ann" \
                 | awk -F/ '{print $3"-"$1"-"$2}' | sort | tail -1)
    fi
    printf "ファイル更新 %s / annotation最新 %s / %6s行  %s\n" \
        "$(date -r "$f" +%Y-%m-%d)" "$latest" "$(wc -l < "$f")" "$d"
done

echo
echo "annotation の最新日付より前に取得することはあり得ないので、"
echo "これが取得日の下限になる。sources.js の acquired に反映すること。"
