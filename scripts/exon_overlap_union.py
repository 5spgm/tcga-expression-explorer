#!/usr/bin/env python3
"""Exonic overlap of every gene with the union of all other genes' exons.

The GDC counters treat reads as strand-agnostic, so the quantity that decides
whether a read is ambiguous is the union of the exons of *all* other genes,
irrespective of strand. This script computes, for every gene in a GENCODE
annotation:

  exon_cov_union  - % of the gene's exonic bases that are exonic in >=1 other gene
  exon_cov_same   - the same restricted to partners on the same strand
  exon_cov_opp    - the same restricted to partners on the opposite strand
  top_partner     - the single partner covering the most of the gene

Usage: exon_overlap_union.py <gencode.gtf[.gz]> <out.tsv>
"""
import gzip
import sys
from collections import defaultdict

BIN = 10_000


def open_any(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def merge(iv):
    iv = sorted(iv)
    out = []
    for s, e in iv:
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def parse(path):
    genes, exons = {}, defaultdict(list)
    with open_any(path) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t")
            if f[2] not in ("gene", "exon"):
                continue
            a = f[8]
            gid = a.split('gene_id "', 1)[1].split('"', 1)[0].split(".")[0]
            if f[2] == "gene":
                nm = a.split('gene_name "', 1)[1].split('"', 1)[0] if 'gene_name "' in a else ""
                gt = a.split('gene_type "', 1)[1].split('"', 1)[0] if 'gene_type "' in a else ""
                genes[gid] = {"chrom": f[0], "strand": f[6], "symbol": nm, "biotype": gt}
            else:
                exons[gid].append((int(f[3]), int(f[4])))
    for gid, xs in exons.items():
        if gid in genes:
            genes[gid]["exons"] = merge(xs)
    return {g: v for g, v in genes.items() if "exons" in v}


def covered(own, others):
    """bases of `own` (merged intervals) covered by `others` (merged intervals)"""
    if not others:
        return 0
    others = merge([list(x) for x in others])
    total, j = 0, 0
    for s, e in own:
        while j < len(others) and others[j][1] < s:
            j += 1
        k = j
        while k < len(others) and others[k][0] <= e:
            total += max(0, min(e, others[k][1]) - max(s, others[k][0]) + 1)
            k += 1
    return total


def main(gtf, out):
    genes = parse(gtf)
    index = defaultdict(list)                      # (chrom, bin) -> [gid]
    for gid, v in genes.items():
        for s, e in v["exons"]:
            for b in range(s // BIN, e // BIN + 1):
                index[(v["chrom"], b)].append(gid)
    rows = []
    for gid, v in genes.items():
        own = v["exons"]
        own_len = sum(e - s + 1 for s, e in own)
        cand = set()
        for s, e in own:
            for b in range(s // BIN, e // BIN + 1):
                cand.update(index[(v["chrom"], b)])
        cand.discard(gid)
        per_partner, same, opp, allx = {}, [], [], []
        for p in cand:
            w = genes[p]
            if w["chrom"] != v["chrom"]:
                continue
            c = covered(own, w["exons"])
            if c <= 0:
                continue
            per_partner[p] = c
            allx += w["exons"]
            (same if w["strand"] == v["strand"] else opp).append(w)
        cov_all = covered(own, allx) if allx else 0
        cov_same = covered(own, [x for w in same for x in w["exons"]]) if same else 0
        cov_opp = covered(own, [x for w in opp for x in w["exons"]]) if opp else 0
        top = max(per_partner, key=per_partner.get) if per_partner else ""
        rows.append((gid, v["symbol"], v["biotype"], v["strand"], own_len,
                     round(100 * cov_all / own_len, 3),
                     round(100 * cov_same / own_len, 3),
                     round(100 * cov_opp / own_len, 3),
                     len(per_partner),
                     genes[top]["symbol"] if top else "",
                     genes[top]["strand"] if top else "",
                     round(100 * per_partner[top] / own_len, 3) if top else 0.0))
    with open(out, "w") as fh:
        fh.write("gene_id\tsymbol\tbiotype\tstrand\texonic_bp\texon_cov_union\texon_cov_same\t"
                 "exon_cov_opp\tn_partners\ttop_partner\ttop_partner_strand\ttop_partner_frac\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print(out, len(rows), "genes", flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
