# TCGA Expression Explorer

A static web tool that shows the expression of one gene, in one TCGA cohort, under
**all three generations of TCGA gene-level RNA-seq quantification** side by side,
on the same specimens, optionally split by published molecular-subtype schemes.

- **Tool:** https://5spgm.github.io/tcga-expression-explorer/
- **Data:** Zenodo [doi:10.5281/zenodo.22670234](https://doi.org/10.5281/zenodo.22670234) (CC BY 4.0)
- **Licence:** MIT for the code in this repository (see `LICENSE`); CC BY 4.0 for the deposited data

| Generation | Period | Pipeline | Annotation | Unit shown |
|---|---|---|---|---|
| Legacy DCC | 2010–2016 | MapSplice + RSEM (hg19) | UCSC hg19 gene model | `normalized_count` |
| GDC v1 (HTSeq era, up to DR-31.0) | 2016–2022 | STAR 2-pass + `htseq-count -m intersection-nonempty -s no` | GENCODE v22 | FPKM-UQ (FPKM for BRCA) |
| GDC v2 (STAR-Counts era, DR-32.0 onwards) | 2022– | STAR `--quantMode GeneCounts`, unstranded | GENCODE v36 | TPM |

## Repository layout

```
docs/                     the web application served by GitHub Pages (plain JS + Plotly)
scripts/preprocess_tcga.py        build the per-gene JSON documents the site fetches
scripts/make_star_matrix.py       assemble a cohort TPM matrix from GDC STAR-Counts files
scripts/make_ffpe_exclude_list.py build the FFPE exclusion list from a GDC sample sheet
scripts/make_*_subtype_table.py   build the subtype tables from the marker-paper supplements
scripts/compare_generations.py    the cross-generation comparison reported in the manuscript
scripts/exon_overlap.py           exonic-overlap fractions from a GENCODE GTF
scripts/make_zenodo_package.py    assemble the deposited archive and its checksums
test/                     a reduced, self-contained test case with expected outputs
PIPELINE.md               how each input file was obtained, in order
```

## Reproducing the analysis in the manuscript

Two commands. Both run on a laptop; neither needs the raw sequencing data.

```bash
pip install -r requirements.txt

# 1. exonic-overlap fractions for every gene in the comparison (GENCODE v36 GTF)
python scripts/exon_overlap.py            # writes exon_overlap_v36_all.tsv

# 2. the three-generation comparison, on the specimens the tool displays
python scripts/compare_generations.py \
    --matrix-dir      <unpacked Zenodo deposit>/matrix \
    --exclude-dir     reference/exclude \
    --nested-stranded reference/overlap.gene.stranded.tsv \
    --nested-strandless reference/overlap.gene.strandless.tsv \
    --exon-overlap    exon_overlap_v36_all.tsv \
    --out-dir         comparison_out
```

`comparison_out/summary.json` holds the per-cohort counts, percentages, odds ratios and
*P* values quoted in the manuscript; `comparison_out/genes_<COHORT>.tsv` holds the
per-gene detection fractions and medians behind them.

Inputs that are not in this repository:

- the expression matrices — download the Zenodo deposit above
- `gencode.v36.annotation.gtf.gz` — https://www.gencodegenes.org/human/release_36.html
- `overlap.gene.stranded.tsv`, `overlap.gene.strandless.tsv` — distributed by the GDC with the
  [mRNA Analysis Pipeline documentation](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/)

## Test case with expected output

`test/` contains a reduced version of one cohort (COAD): the three generation matrices cut
down to a gene subset and a specimen subset, the exclusion list, and the expected outputs.

```bash
bash test/run_test.sh
```

The script runs `scripts/compare_generations.py` on `test/input/` and compares the result with
`test/expected/`, printing `TEST PASSED` when every value matches. It takes a few seconds.
The test includes *TREX1* (nested inside *ATRIP*; zero in both GDC generations), *ATRIP*
(the host gene, unaffected) and a set of non-nested controls, so a failure localises
immediately to either the specimen filter, the identifier mapping or the detection rule.

## Rebuilding the site data

`PIPELINE.md` documents the order; `scripts/run_production.sh` and the per-cohort
`scripts/run_<cohort>.sh` wrappers drive it. Specimen selection is identical in the site and
in the comparison: sample-type codes 01/03/05/09 are tumour, codes ≥ 10 are normal, codes
02/04/06/07/08/40 (recurrent and metastatic) are excluded by default, duplicate vials of the
same patient and sample type are collapsed to the first vial, and the per-cohort exclusion
lists remove FFPE specimens (and, for PAAD, non-PDAC histology and pseudonormal specimens).

## Exonic overlap (primary classifier)

`scripts/exon_overlap_union.py <gencode.gtf> <out.tsv>` computes, for every gene in a
GENCODE annotation, the share of its exonic bases that are also exonic in **at least one
other gene** — the union over all partners and both strands, which is the quantity that
decides ambiguity for the GDC's unstranded counters. It also reports the same share
restricted to same-strand and to opposite-strand partners, and the single partner covering
the most of the gene. Run it once per annotation:

```
python scripts/exon_overlap_union.py gencode.v36.annotation.gtf exon_cov_union_v36.tsv
python scripts/exon_overlap_union.py gencode.v22.annotation.gtf exon_cov_union_v22.tsv
```

`scripts/compare_generations.py --exon-overlap exon_cov_union_v36.tsv` then classifies each
generation by the annotation it was counted against. (`scripts/exon_overlap.py` is the earlier,
per-partner version, kept because the first submission used it.)

## Direct test of the counting rule

`scripts/run_validation.sh` downloads nothing itself: give it two paired FASTQ files in
`validation/fastq/` (the manuscript uses ENCODE experiments ENCSR985KAT and ENCSR999ZCI),
a STAR index and GENCODE v36/v22 GTFs, and it aligns each library once and counts the same
alignment six ways — STAR `--quantMode GeneCounts`, `htseq-count -m intersection-nonempty`
against v36 and v22, `htseq-count -m union`, `featureCounts` default, and
`featureCounts -O --fraction`. `scripts/extract_validation.py` collects *TREX1*, *ATRIP*,
*SHISA5* and *ACTB* from the six outputs together with each counter's ambiguous-read tally
and writes `validation_counts.tsv` (Supplementary Table S5 of the manuscript).

Requires STAR 2.7.11b, HTSeq 2.0.5 and Subread 2.0.6 in addition to `requirements.txt`;
both are available from bioconda.

## Citation

If you use the tool or the data, please cite the manuscript (in preparation) and the Zenodo
deposit. A machine-readable `CITATION.cff` will be added when the article has a DOI.
