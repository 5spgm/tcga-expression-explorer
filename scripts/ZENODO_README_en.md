# TCGA Expression Explorer — data underlying the published figures

Curated RNA-seq matrices and sample annotation used to compare gene expression
across three generations of TCGA quantification pipelines. This is the exact
input to the **TCGA Expression Explorer** web tool; every figure the tool
produces can be regenerated from these files.

- Tool: https://5spgm.github.io/tcga-expression-explorer/
- Code: https://github.com/5spgm/tcga-expression-explorer

---

## Why deposit these data

TCGA RNA-seq has been requantified three times. **The intermediate (HTSeq)
generation was replaced by STAR-Counts at Data Release 32 in 2022, and the same
values can no longer be retrieved from the current GDC.**

| Generation | Pipeline | Value | Gene rows | Obtained from |
|---|---|---|---|---|
| Old | UNC IlluminaHiSeq RNASeqV2 | normalized count | 20,531 | Former TCGA Data Coordinating Center (http://cancergenome.nih.gov/, closed) |
| Mid | GDC HTSeq (DR15–31) | FPKM / FPKM-UQ | 60,483 | The GDC at the time. **Absent from the current GDC** |
| New | GDC STAR-Counts (DR32+, 2022–) | TPM | 60,660 | GDC (current) |

The GDC replaces earlier quantifications at each data release, so the HTSeq-era
matrices survive only where individual users downloaded them between DR15 and
DR31. Preserving the ability to audit how expression values changed across
these transitions is the main motivation for this deposit.

There is one exception for the intermediate generation: the breast cancer
intermediate matrix in this package was downloaded from the UCSC Xena GDC hub on
**2024-04-13**, as `TCGA-BRCA.htseq_fpkm.tsv.gz` (the date is the timestamp of
the local .gz; the gzip header still carries the distribution-side file time of
2019-07-18 22:22:54 UTC, and the SHA-256 of the extracted tsv matches the stream decompressed
from that .gz). What Xena distributes is log2(FPKM+1)
rather than raw FPKM (see Contents below). For every other cohort, a copy
downloaded at the time is the only route to these values.

The hub catalogue has since moved to the STAR era (gencode v36) and the
**cohort-level HTSeq FPKM datasets (`TCGA-BRCA.htseq_fpkm`) have been removed.**
As of 2026-09-08 the only HTSeq FPKM left on the hub is the pan-cancer
`GDC-PANCAN.htseq_fpkm` (gencode v22). All 1,217 breast cancer samples in this
package are present in it, and the values agree to the precision the hub's API
returns. The same values therefore remain obtainable by this route:

    hub      https://gdc.xenahubs.net
    dataset  GDC-PANCAN.htseq_fpkm.tsv   (gencode v22, log2(FPKM+1))
    page     https://xenabrowser.net/datapages/?dataset=GDC-PANCAN.htseq_fpkm.tsv
             &host=https%3A%2F%2Fgdc.xenahubs.net
    method   take the columns of the breast cancer samples out of the
             pan-cancer matrix (sample list: sample_annotation/BRCA_samples.tsv)

The old generation was obtained **in its entirety from the former TCGA portal
(http://cancergenome.nih.gov/)**, which has since been closed, with TCGA data
migrated to the NIH GDC. Equivalent values are redistributed by UCSC Xena and
others, but with gene symbols alone as identifiers, unlike the
"SYMBOL|EntrezID" retained here. The new generation is obtainable from the
current GDC. Both are included so that the exact input behind the published
figures is preserved.

---

## Contents

```
matrix/<cancer>/<cancer>_<generation>_<value_type>.tsv.gz
    Expression matrices, kept in the form in which they were obtained.
    They are deliberately NOT merged into a single table: gene identifier
    systems and sample sets differ between generations, so a merged matrix
    would be mostly missing values, and the identifier conversion — our
    interpretation rather than source data — would become unverifiable.

      old  first column "SYMBOL|EntrezID", 20,531 rows (hg19, pre-GENCODE)
      mid  first column Ensembl gene ID, 60,483 rows (GENCODE v22)
      new  gene_id / gene_name / gene_type columns + samples,
           60,660 rows (GENCODE v36)

    The row counts are identical across all five cohorts, so a file can be
    assigned to its generation on its own.

    The breast cancer intermediate matrix alone holds log2(FPKM+1). It is
    placed here exactly as the UCSC Xena GDC hub distributes it, without
    back-transformation, and the filename says so:
    BRCA_mid_FPKM_log2p1.tsv.gz. Apply 2^x - 1 to recover raw FPKM. The
    intermediate matrices of the other cohorts are untransformed FPKM-UQ.

    This single file is the .gz obtained from UCSC Xena, placed here
    byte-for-byte with no recompression and no rewriting (SHA-256
    b5856c84…). Two consequences distinguish it from the other 14 files:
      - its first column is named Ensembl_ID, not gene_id
      - its gzip header still carries the distribution-side filename
        TCGA-BRCA.htseq_fpkm.tsv and file time 2019-07-18 22:22:54 UTC
        (gzip stores UTC seconds in the header and gzip -l renders them in the
        local timezone, so it prints 2019-07-19 07:22:54 in JST)
        (inspect with gzip -l or zcat -f)
    The figures on the public site were produced from an xlsx conversion of
    this same .gz. Across the first 100 genes x all 1,217 samples the two
    agree to within 3.2e-14 (double-precision rounding), well below any
    displayed digit.

sample_annotation/<cancer>_samples.tsv
    One row per sample: which generations contain it, which group it belongs
    to, and which subtype assignments apply. Barcode granularity differs
    between generations (the old generation omits the vial letter), so rows
    are keyed on patient ID plus sample type code.

exclude/
    Samples removed from the analysis, with the basis for removal.

provenance/<cancer>/<generation>/
    The GDC sample sheet and download manifest as they were at the time of
    download, and for the old generation the TCGA DCC bundle
    (file_manifest.txt, FILE_SAMPLE_MAP.txt, README_DCC.txt, and for stomach
    and pancreatic cancer also file_annotations.txt).
    file_annotations.txt holds the annotations the DCC attached, which cannot
    be consulted in the same form on the current GDC. Its content differs
    between the two cohorts:
      stomach     aliquot-level Item flagged DNU records
                  (1,212 rows = 11 aliquots) — flagged by the DCC as
                  not to be used for analysis
      pancreatic  patient-level notes of study-protocol deviations (10 rows) The GDC replaces its distribution at every data release, so
    these files are the only record of which files existed and were retrieved
    on that date — in particular for the HTSeq-era matrices that can no longer
    be obtained. Download dates are recorded in inventory.json.

    14 of the 15 layers (5 cohorts x 3 generations) are covered. The only
    one absent is the breast cancer intermediate layer: it was not downloaded
    from the GDC (it comes from the UCSC Xena GDC hub), so no GDC sample
    sheet exists to deposit.

inventory.json      list of deposited matrices (dimensions, source filenames)
checksums.sha256    SHA-256 for every file
```

### Columns of sample_annotation

| Column | Meaning |
|---|---|
| `patient_id` | TCGA-XX-XXXX |
| `sample_type_code` | 01 primary, 02 recurrent, 06 metastatic, 11 solid tissue normal, etc. **Zero-padded string** |
| `group` | Tumor (primary) / TumorExtra (recurrent or metastatic) / Normal |
| `generations` | generations in which the sample is present |
| `barcode_old` / `_mid` / `_new` | the actual barcode in each generation |
| `subtype_*` | subtype assignments (see references below) |

Read `sample_type_code` as a string:

```python
pd.read_csv("PAAD_samples.tsv", sep="\t", dtype={"sample_type_code": str})
```

---

## Sample selection

Samples were removed before any figure was produced. Every decision can be
checked against `exclude/` and `sample_annotation/`.

- **Recurrent and metastatic samples are separated from primary tumours.**
  Pooling them by default risks silent misinterpretation; the `group` column
  distinguishes them.
- **Samples flagged by the DCC are not excluded.** The 11 DNU-flagged
  aliquots in stomach cancer are present in all three generations (11 samples
  in the old matrix, 9 in the intermediate, 11 in the new), and the 10
  patients carrying study-protocol deviation notes in pancreatic cancer are
  present in all three as well. The published figures were produced with
  these samples included. To exclude them, use
  `provenance/STAD/old/file_annotations.txt` and
  `provenance/PAAD/old/file_annotations.txt`.

- **FFPE-derived samples are excluded.** The vial letter (01A / 01B) does not
  indicate FFPE status, so exclusion is based on the GDC metadata fields
  `is_ffpe` / `Preservation Method`. In breast cancer 13 samples qualified,
  all listed in `exclude/`. **One of them, `TCGA-PL-A8LZ-01A`, carries the
  01A vial letter** — precisely why the vial letter cannot serve as a proxy.
  How many are actually dropped depends on the generation: all 13 in the
  intermediate matrix (1,217 samples), 2 in the new (1,208) and none in the
  old (1,154). The remaining 11 are 01B / 01C vials, and the new and old
  matrices carry only the 01A vial of those samples, so they were never
  present to begin with.
- **Duplicate vials of the same sample are collapsed to one**, otherwise the
  same patient is counted twice. The lowest vial letter is retained.
- **In pancreatic cancer, 28 samples are excluded**: 19 determined not to be
  pancreatic ductal adenocarcinoma, and 9 with very low neoplastic cellularity
  (Raphael et al. 2017, Table S1).

---

## Sources

Expression data were obtained from:

- **Former TCGA Data Coordinating Center** (http://cancergenome.nih.gov/) —
  the old generation (RNASeqV2 normalized count) for all five cohorts. This
  portal has been closed and TCGA data migrated to the NIH GDC.
- **NCI Genomic Data Commons** (https://portal.gdc.cancer.gov/) — the
  intermediate generation (HTSeq FPKM / FPKM-UQ, except breast cancer) and the
  new generation (STAR-Counts TPM).
- **UCSC Xena** (https://xenabrowser.net/) — breast cancer intermediate only,
  via the HTSeq FPKM (log2(FPKM+1)) redistributed by its GDC hub.

**Use of TCGA data is subject to the terms of the respective data providers.**

Subtype assignments (`subtype_*` columns) have two kinds of origin. **Please
cite the original publication or data provider in each case.**

Triple-negative status (`subtype_TNBC`) is derived from the ER / PR / HER2
calls in the TCGA BCR biospecimen clinical file
(`nationwidechildrens.org_clinical_patient_brca.txt`), not from the
supplementary material of a paper.

The following derive from published supplementary material.

- Hoadley KA, et al. *Cell* 2018 — pan-cancer iCluster
- Liu Y, et al. *Cancer Cell* 2018 — gastrointestinal adenocarcinoma subtypes, MSI, CIMP, CMS
- Thorsson V, et al. *Immunity* 2018 — immune subtypes C1–C6, published per-cancer subtypes including PAM50
- Raphael BJ, et al. *Cancer Cell* 2017 — pancreatic ductal adenocarcinoma (Moffitt / Bailey / Collisson / purity)

---

## Reproducing the figures

Preprocessing scripts are in `scripts/` of the GitHub repository. Using
`matrix/` from this deposit as input reproduces the JSON served by the site.

```bash
python3 scripts/preprocess_tcga.py \
    --cancer-type PAAD \
    --new-tpm       matrix/PAAD/PAAD_new_TPM.tsv.gz \
    --mid-fpkmuq    matrix/PAAD/PAAD_mid_FPKM_UQ.tsv.gz \
    --old-normcount matrix/PAAD/PAAD_old_normalized_count.tsv.gz \
    --exclude-samples exclude/PAAD_exclude_all.txt \
    --out-dir ./data
```

For breast cancer the intermediate matrix is on a different scale, so pass
**both** `--mid-fpkm` and `--mid-log2p1`. Omitting `--mid-log2p1` would mean
reading log2 values as FPKM; the script stops with an error in that case, so
it cannot silently produce wrong figures.

```bash
python3 scripts/preprocess_tcga.py \
    --cancer-type BRCA \
    --new-tpm       matrix/BRCA/BRCA_new_TPM.tsv.gz \
    --mid-fpkm      matrix/BRCA/BRCA_mid_FPKM_log2p1.tsv.gz \
    --mid-log2p1 \
    --old-normcount matrix/BRCA/BRCA_old_normalized_count.tsv.gz \
    --exclude-samples exclude/BRCA_ffpe_exclude_brca.txt \
    --out-dir ./data
```

See `PIPELINE.md` in the repository for the full procedure.

## Verifying integrity

```bash
sha256sum -c checksums.sha256
```

## License

The organisation and annotation of this package are released under CC-BY 4.0.
The underlying TCGA data and the supplementary material of the cited papers
remain subject to the terms of their respective providers.
