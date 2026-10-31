#!/usr/bin/env bash
# Reduced end-to-end test for scripts/compare_generations.py.
# Runs the comparison on test/input (one cohort, 294 genes, 20 specimens) and
# compares the result with test/expected. Takes a few seconds.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
OUT="$(mktemp -d)"
trap 'rm -rf "$OUT"' EXIT

python "$ROOT/scripts/compare_generations.py" \
    --matrix-dir        "$HERE/input" \
    --exclude-dir       "$HERE/input/exclude" \
    --nested-stranded   "$HERE/reference/overlap.gene.stranded.tsv" \
    --nested-strandless "$HERE/reference/overlap.gene.strandless.tsv" \
    --exon-overlap      "$HERE/reference/exon_cov_union_v36_subset.tsv" \
    --cohorts           COAD \
    --out-dir           "$OUT" > "$OUT/stdout.txt"

fail=0
for f in genes_COAD.tsv summary.json; do
    if diff -q "$HERE/expected/$f" "$OUT/$f" > /dev/null; then
        echo "  [ok]     $f"
    else
        echo "  [DIFFER] $f"
        diff "$HERE/expected/$f" "$OUT/$f" | head -20
        fail=1
    fi
done
if [ "$fail" -eq 0 ]; then
    echo "TEST PASSED"
else
    echo "TEST FAILED"; exit 1
fi
