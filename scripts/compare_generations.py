#!/usr/bin/env python3
"""
compare_generations.py
======================
Compare gene detection across the three generations of TCGA gene-level RNA-seq
(legacy DCC / GDC HTSeq era / GDC STAR-Counts era) on the *same specimens*,
using exactly the specimen selection the TCGA Expression Explorer displays.

Specimen selection (identical to scripts/preprocess_tcga.py):
  * sample-type code 01/03/05/09 -> Tumour, codes >= 10 -> Normal,
    codes 02/04/06/07/08/40 (recurrent / metastatic) -> excluded
  * explicit exclusion list per cohort (FFPE from the GDC sample sheet; for PAAD
    additionally non-PDAC histology and pseudonormal specimens)
  * duplicate vials of the same patient x sample-type collapsed to the first vial
  * the comparison is then restricted to specimens present in all three generations

Gene classification:
  * GDC "Overlapped Genes" lists (stranded / strandless) -> nested same strand /
    nested opposite strand / not nested
  * exonic overlap fraction computed from the GENCODE v36 GTF (exon_overlap.py)

Outputs (per cohort and pooled): gene-level tables, a summary JSON and the
Fisher tests reported in the manuscript.

    python compare_generations.py --matrix-dir <zenodo>/matrix --out-dir out
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

COHORTS = ["BRCA", "COAD", "LUAD", "PAAD", "STAD"]
MID_UNIT = {"BRCA": "FPKM_log2p1", "COAD": "FPKM_UQ", "LUAD": "FPKM_UQ",
            "PAAD": "FPKM_UQ", "STAD": "FPKM_UQ"}
TUMOUR_CODES = {"01", "03", "05", "09"}
EXTRA_CODES = {"02", "04", "06", "07", "08", "40"}
BARCODE_RE = re.compile(r"^(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})-(\d{2})([A-Z]?)")
DET_LEGACY = 0.5      # detected in legacy: non-zero in more than half the specimens
LOST_GDC = 0.1        # lost: non-zero in at most 10% in both GDC generations


def classify(barcode: str):
    """-> (patient, type_code, vial) or None when the column is not a specimen."""
    m = BARCODE_RE.match(str(barcode).strip())
    if not m:
        return None
    return m.group(1), m.group(2), m.group(3)


def select_specimens(columns, excluded: set[str]):
    """Apply the tool's specimen filter. Returns {(patient, code): column}."""
    keep: dict[tuple[str, str], list[tuple[str, str]]] = defaultdict(list)
    for col in columns:
        parsed = classify(col)
        if parsed is None:
            continue
        patient, code, vial = parsed
        if code in EXTRA_CODES:
            continue
        if code not in TUMOUR_CODES and int(code) < 10:
            continue
        if col.strip() in excluded or patient in excluded:
            continue
        keep[(patient, code)].append((vial, col))
    return {k: sorted(v)[0][1] for k, v in keep.items()}


def read_matrix(path: Path, id_kind: str):
    """Read a deposited matrix. Returns (DataFrame indexed by key, gene_name map)."""
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        first_row = fh.readline().rstrip("\n").split("\t")
    if len(first_row) == len(header) + 1:
        # no corner label (R write.table default): the first data column is the gene id
        df = pd.read_csv(path, sep="\t", index_col=0, low_memory=False)
        df.index = df.index.astype(str)
        if id_kind == "mid":
            return df.rename(index=lambda s: re.sub(r"\.\d+$", "", s)), None
        if id_kind == "old":
            return df.rename(index=lambda s: s.split("|")[0]), None
        raise ValueError(f"unexpected header shift in {path}")
    # keep column 0 even when it has no name (an empty corner label is the id column);
    # drop any other unnamed column (trailing tab in some legacy exports)
    usecols = [0] + [i for i, c in enumerate(header) if i > 0 and c.strip() != ""]
    df = pd.read_csv(path, sep="\t", usecols=usecols, low_memory=False)
    first = df.columns[0]
    names = None
    if id_kind == "new":
        name_col = next((c for c in df.columns[:4] if c.lower() in ("gene_name", "symbol")), None)
        type_col = next((c for c in df.columns[:4] if c.lower() == "gene_type"), None)
        key = df[first].astype(str).str.replace(r"\.\d+(_PAR_Y)?$", "", regex=True)
        names = pd.DataFrame({"ensembl": key, "symbol": df[name_col] if name_col else None,
                              "gene_type": df[type_col] if type_col else None})
        drop = [c for c in (first, name_col, type_col) if c]
        df = df.drop(columns=drop)
        df.index = key
    elif id_kind == "mid":
        key = df[first].astype(str).str.replace(r"\.\d+$", "", regex=True)
        df = df.drop(columns=[first]); df.index = key
    else:  # legacy: "SYMBOL|EntrezID"
        key = df[first].astype(str).str.split("|").str[0]
        df = df.drop(columns=[first]); df.index = key
    return df, names


def cohort_table(mat_dir: Path, cohort: str, excluded: set[str], nested: dict, exon_cov: dict):
    old, _ = read_matrix(mat_dir / cohort / f"{cohort}_old_normalized_count.tsv.gz", "old")
    mid, _ = read_matrix(mat_dir / cohort / f"{cohort}_mid_{MID_UNIT[cohort]}.tsv.gz", "mid")
    new, names = read_matrix(mat_dir / cohort / f"{cohort}_new_TPM.tsv.gz", "new")
    if MID_UNIT[cohort] == "FPKM_log2p1":
        mid = np.power(2.0, mid) - 1.0

    sel = [select_specimens(d.columns, excluded) for d in (old, mid, new)]
    common = sorted(set(sel[0]) & set(sel[1]) & set(sel[2]))
    cols = [[s[k] for k in common] for s in sel]
    codes = sorted({k[1] for k in common})

    # legacy rows are gene symbols; map the GDC generations onto symbols via gene_name
    ens2sym = dict(zip(names["ensembl"], names["symbol"].astype(str)))
    sym_counts = pd.Series(list(ens2sym.values())).value_counts()
    unique_syms = set(sym_counts[sym_counts == 1].index)

    new_sym = new.loc[:, cols[2]].copy()
    new_sym["symbol"] = [ens2sym.get(i) for i in new_sym.index]
    new_sym = new_sym[new_sym["symbol"].isin(unique_syms)].set_index("symbol")
    new_sym = new_sym[~new_sym.index.duplicated(keep="first")]
    mid_sym = mid.loc[:, cols[1]].copy()
    mid_sym["symbol"] = [ens2sym.get(i) for i in mid_sym.index]
    mid_sym = mid_sym[mid_sym["symbol"].isin(unique_syms)].set_index("symbol")
    mid_sym = mid_sym[~mid_sym.index.duplicated(keep="first")]
    old_sym = old.loc[:, cols[0]]
    old_sym = old_sym[~old_sym.index.duplicated(keep="first")]

    shared = sorted(set(old_sym.index) & set(mid_sym.index) & set(new_sym.index))
    o = old_sym.loc[shared].to_numpy(dtype=float)
    m = mid_sym.loc[shared].to_numpy(dtype=float)
    n = new_sym.loc[shared].to_numpy(dtype=float)

    def frac_nonzero(a):
        with np.errstate(invalid="ignore"):
            return np.nanmean(np.nan_to_num(a, nan=0.0) > 0, axis=1)

    tab = pd.DataFrame({
        "symbol": shared,
        "det_old": frac_nonzero(o), "det_mid": frac_nonzero(m), "det_new": frac_nonzero(n),
        "med_old": np.nanmedian(o, axis=1), "med_mid": np.nanmedian(m, axis=1),
        "med_new": np.nanmedian(n, axis=1),
    })
    sym2ens = {v: k for k, v in ens2sym.items() if v in unique_syms}
    tab["gene_id"] = [sym2ens.get(s, "") for s in tab["symbol"]]
    tab["group"] = [nested.get(s, "not nested") for s in tab["symbol"]]
    tab["exon_cov"] = [exon_cov.get(g, 0.0) for g in tab["gene_id"]]
    tab["legacy_detected"] = tab["det_old"] > DET_LEGACY
    tab["lost_in_gdc"] = tab["legacy_detected"] & (tab["det_mid"] <= LOST_GDC) & (tab["det_new"] <= LOST_GDC)
    return tab, len(common), len(shared), codes


def enrichment(tab: pd.DataFrame, key: str, positive) -> dict:
    det = tab[tab["legacy_detected"]]
    pos = det[positive(det)]
    neg = det[~positive(det)]
    a, b = int(pos["lost_in_gdc"].sum()), int(neg["lost_in_gdc"].sum())
    orr, pv = fisher_exact([[a, len(pos) - a], [b, len(neg) - b]])
    return {"criterion": key, "n_pos": len(pos), "lost_pos": a,
            "pct_pos": 100 * a / len(pos) if len(pos) else None,
            "n_neg": len(neg), "lost_neg": b,
            "pct_neg": 100 * b / len(neg) if len(neg) else None,
            "odds_ratio": float(orr), "p_value": float(pv)}


def threshold_sensitivity(tab: pd.DataFrame) -> dict:
    """Repeat the >=90% enrichment under other definitions of detected / lost."""
    out = {}
    for det_thr, lost_thr in [(0.5, 0.1), (0.8, 0.02), (0.9, 0.0), (0.3, 0.2)]:
        t = tab.copy()
        t["legacy_detected"] = t["det_old"] > det_thr
        t["lost_in_gdc"] = t["legacy_detected"] & (t["det_mid"] <= lost_thr) & (t["det_new"] <= lost_thr)
        e = enrichment(t, f"detected >{det_thr:g}, lost <={lost_thr:g}", lambda d: d["exon_cov"] >= 90)
        out[f"det>{det_thr:g}_lost<={lost_thr:g}"] = e
    return out


def lost_vs_kept(tab: pd.DataFrame, exon_meta: pd.DataFrame) -> dict:
    """Compare the genes with >=90% exonic overlap that are lost with those that are not."""
    det = tab[tab["legacy_detected"]]
    hi = det[det["exon_cov"] >= 90]
    res = {}
    for label, d in [("lost", hi[hi["lost_in_gdc"]]), ("kept", hi[~hi["lost_in_gdc"]])]:
        meta = exon_meta.reindex(d["gene_id"])
        res[label] = {
            "n": int(len(d)),
            "median_legacy_median": float(d["med_old"].median()) if len(d) else None,
            "median_legacy_detection": float(d["det_old"].median()) if len(d) else None,
            "median_exon_cov": float(d["exon_cov"].median()) if len(d) else None,
            "pct_on_gdc_nested_list": float(100 * (d["group"] != "not nested").mean()) if len(d) else None,
            "median_exonic_bp": float(meta["exonic_bp"].median()) if len(d) else None,
            "median_n_partners": float(meta["n_partners"].median()) if len(d) else None,
            "biotypes": meta["biotype"].value_counts().head(6).to_dict() if len(d) else {},
        }
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix-dir", required=True)
    ap.add_argument("--exclude-dir", required=True, help="directory holding the per-cohort exclusion lists")
    ap.add_argument("--nested-stranded", required=True, help="GDC overlap.gene.stranded.tsv")
    ap.add_argument("--nested-strandless", required=True, help="GDC overlap.gene.strandless.tsv")
    ap.add_argument("--exon-overlap", required=True, help="exon_overlap_v36_all.tsv from exon_overlap.py")
    ap.add_argument("--cohorts", default=",".join(COHORTS),
                    help="comma-separated cohorts to process (default: all five)")
    ap.add_argument("--out-dir", default="comparison_out")
    args = ap.parse_args()

    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    st = pd.read_csv(args.nested_stranded, sep="\t")
    un = pd.read_csv(args.nested_strandless, sep="\t")
    nested = {s: "nested same strand" for s in st["gene_name"].astype(str)}
    for s in un["gene_name"].astype(str):
        nested.setdefault(s, "nested opposite strand")
    ex = pd.read_csv(args.exon_overlap, sep="\t")
    if "exon_cov_union" in ex.columns:
        # strand-agnostic union of every other gene's exons: the quantity that
        # decides ambiguity for an unstranded counter
        ex["cov"] = ex["exon_cov_union"].fillna(0)
    else:
        ex["cov"] = ex[["exonic_frac_same_strand_max", "exonic_frac_opp_strand_max"]].max(axis=1).fillna(0)
    # key the coverage on the Ensembl gene id: GENCODE v36 carries 1,189 duplicated
    # gene symbols, and the comparison resolves symbols through the GDC v2 matrix
    ex["gene_id"] = ex["gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)
    exon_cov = dict(zip(ex["gene_id"], ex["cov"]))
    exon_meta = ex.drop_duplicates("gene_id").set_index("gene_id")

    summary = {}
    for cohort in [c.strip().upper() for c in args.cohorts.split(",") if c.strip()]:
        f = Path(args.exclude_dir) / f"exclude_{cohort.lower()}.txt"
        excluded = {l.strip() for l in f.read_text().splitlines() if l.strip()} if f.exists() else set()
        tab, n_spec, n_genes, codes = cohort_table(Path(args.matrix_dir), cohort, excluded, nested, exon_cov)
        tab.to_csv(out / f"genes_{cohort}.tsv", sep="\t", index=False)
        det = tab[tab["legacy_detected"]]
        per_group = {g: {"n": int((det["group"] == g).sum()),
                         "lost": int(det.loc[det["group"] == g, "lost_in_gdc"].sum())}
                     for g in ["nested same strand", "nested opposite strand", "not nested"]}
        for g in per_group:
            per_group[g]["pct"] = (100 * per_group[g]["lost"] / per_group[g]["n"]) if per_group[g]["n"] else None
        strict = det[det["med_old"] > 1]
        summary[cohort] = {
            "n_specimens_compared": n_spec, "n_specimens_excluded": len(excluded),
            "sample_type_codes": codes,
            "n_genes_compared": n_genes, "n_legacy_detected": int(len(det)),
            "by_group": per_group,
            "enrichment_nested": enrichment(tab, "nested (either strand) vs not nested",
                                            lambda d: d["group"] != "not nested"),
            "enrichment_exonic90": enrichment(tab, "exonic overlap >= 90% vs < 90%",
                                              lambda d: d["exon_cov"] >= 90),
            "exonic_bins": {lab: {"n": int(((det["exon_cov"] >= lo) & (det["exon_cov"] < hi)).sum()),
                                  "lost": int(det.loc[(det["exon_cov"] >= lo) & (det["exon_cov"] < hi),
                                                      "lost_in_gdc"].sum())}
                            for lo, hi, lab in [(0, 50, "<50%"), (50, 90, "50-90%"), (90, 100.01, ">=90%")]},
            "sensitivity_legacy_median_gt1": {
                "n_legacy_detected": int(len(strict)),
                "enrichment_exonic90": enrichment(tab[tab["med_old"] > 1], "exonic >= 90%, legacy median > 1",
                                                  lambda d: d["exon_cov"] >= 90)},
            "sensitivity_detection_thresholds": threshold_sensitivity(tab),
            "exonic90_lost_vs_kept": lost_vs_kept(tab, exon_meta),
        }
        e = summary[cohort]["enrichment_exonic90"]
        if e["pct_pos"] is None or e["pct_neg"] is None:
            print(f"{cohort}: {n_spec} specimens, {n_genes} genes, "
                  f"no gene in one of the two exonic-overlap classes", flush=True)
        else:
            print(f"{cohort}: {n_spec} specimens, {n_genes} genes, "
                  f"{e['pct_pos']:.1f}% vs {e['pct_neg']:.3f}% lost "
                  f"(OR {e['odds_ratio']:,.0f}, P {e['p_value']:.1e})", flush=True)

    (out / "summary.json").write_text(json.dumps(summary, indent=1))
    print("wrote", out / "summary.json")


if __name__ == "__main__":
    main()
