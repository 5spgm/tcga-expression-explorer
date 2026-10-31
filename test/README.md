# Test case

`run_test.sh` runs `scripts/compare_generations.py` end to end on a reduced input and
compares the output with the stored expectation.

- `input/COAD/` — the three generation matrices for COAD, cut down to 294 genes and 20
  specimens (12 primary tumours, 8 solid-normal). The three files keep the formatting
  quirks of the originals on purpose: the legacy matrix has an empty corner label, the
  HTSeq-era matrix has a header one field shorter than its rows, and the STAR-Counts
  matrix carries `gene_name` and `gene_type` columns. A reader that silently mishandles
  any of these fails the test.
- `input/exclude/exclude_coad.txt` — empty; COAD has no excluded specimens.
- `reference/` — the two GDC nested-gene lists and the exonic-overlap table restricted
  to the gene subset.
- `expected/` — `genes_COAD.tsv` and `summary.json` as produced by the pinned
  dependencies in `requirements.txt`.

The gene subset is **deliberately enriched** for genes that lose detection (all COAD
losses plus *TREX1*, *ATRIP*, *SHISA5*, two housekeeping genes, 40 nested genes that are
retained and 200 random non-nested genes). The percentages in `expected/summary.json`
are therefore not the percentages reported in the manuscript; they are fixed values that
change only if the code changes. To reproduce the manuscript numbers, run the same script
on the full deposit as shown in the top-level `README.md`.
