import gzip, re, json, sys
import numpy as np, pandas as pd

def parse(path):
    genes, exons = {}, {}
    with gzip.open(path, "rt") as fh:
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
                genes[gid] = (f[0], int(f[3]), int(f[4]), f[6], nm, gt)
            else:
                exons.setdefault(gid, []).append((f[0], int(f[3]), int(f[4]), f[6]))
    return genes, exons

def merge(iv):
    iv = sorted(iv)
    out = []
    for s, e in iv:
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out

def exonic(exons):
    """gene -> {chrom: merged intervals}, strand kept outside"""
    d = {}
    for g, xs in exons.items():
        by = {}
        for c, s, e, st in xs:
            by.setdefault(c, []).append((s, e))
        d[g] = {c: merge(v) for c, v in by.items()}
    return d

def olap_len(a, b):
    """total overlap length between two sorted merged interval lists"""
    i = j = t = 0
    while i < len(a) and j < len(b):
        s = max(a[i][0], b[j][0]); e = min(a[i][1], b[j][1])
        if e >= s:
            t += e - s + 1
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return t

def analyse(gtf, targets, tag):
    genes, exons = parse(gtf)
    name2id = {}
    for g, v in genes.items():
        name2id.setdefault(v[4], g)
    targets = [t if t in genes else name2id.get(t, t) for t in targets]
    ex = exonic(exons)
    # chromosome -> list of (start, end, gid) gene spans for candidate search
    bych = {}
    for g, (c, s, e, st, nm, gt) in genes.items():
        bych.setdefault(c, []).append((s, e, g))
    for c in bych:
        bych[c].sort()
    rows = []
    for g in targets:
        if g not in genes or g not in ex:
            rows.append({"gene_id": g, "found": False, "symbol": g})
            continue
        c, s, e, st, nm, gt = genes[g]
        mine = ex[g].get(c, [])
        mylen = sum(b - a + 1 for a, b in mine)
        same = opp = 0
        partners = []
        for s2, e2, g2 in bych.get(c, []):
            if g2 == g or e2 < s or s2 > e:
                continue
            o = olap_len(mine, ex.get(g2, {}).get(c, []))
            if o <= 0:
                continue
            st2 = genes[g2][3]
            if st2 == st:
                same = max(same, 0) + 0
                same += o
            else:
                opp += o
            partners.append((genes[g2][4], g2, st2, o, round(100 * o / mylen, 1)))
        partners.sort(key=lambda x: -x[3])
        best_same = max([p[4] for p in partners if p[2] == st], default=0.0)
        best_opp = max([p[4] for p in partners if p[2] != st], default=0.0)
        rows.append({"gene_id": g, "found": True, "symbol": nm, "biotype": gt, "strand": st,
                     "exonic_bp": mylen,
                     "exonic_frac_same_strand_max": best_same,
                     "exonic_frac_opp_strand_max": best_opp,
                     "n_partners": len(partners),
                     "top_partner": partners[0][0] if partners else "",
                     "top_partner_strand": partners[0][2] if partners else "",
                     "top_partner_frac": partners[0][4] if partners else 0.0})
    pd.DataFrame(rows).to_csv(f"exon_overlap_{tag}.tsv", sep="\t", index=False)
    print(tag, "done", len(rows), flush=True)

if __name__ == "__main__":
    tg = json.load(open("exon_targets.json"))
    analyse("ref/gencode.v36.annotation.gtf.gz", tg, "v36")
    import os
    if os.path.exists("ref/gencode.v22.annotation.gtf.gz"):
        analyse("ref/gencode.v22.annotation.gtf.gz", tg, "v22")
