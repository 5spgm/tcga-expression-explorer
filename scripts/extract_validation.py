#!/usr/bin/env python3
"""Collect TREX1 / ATRIP / SHISA5 counts from the six counting conditions."""
import json
import os
import re
import sys

import pandas as pd

W = "/media/nighthawk/_mnt_16TB/webtool/validation"
GENES = {"ENSG00000213689": "TREX1", "ENSG00000164053": "ATRIP", "ENSG00000280804": "ATRIP-TREX1 readthrough", "ENSG00000164054": "SHISA5",
         "ENSG00000075624": "ACTB"}
SAMPLES = ["HepG2", "sigmoid_colon"]


def base(gid):
    return gid.split(".")[0]


def read_htseq(path):
    vals, special = {}, {}
    for line in open(path):
        k, v = line.rstrip("\n").split("\t")[:2]
        (special if k.startswith("__") else vals)[base(k)] = float(v)
    return vals, special


def read_fc(path):
    df = pd.read_csv(path, sep="\t", comment="#")
    col = df.columns[-1]
    df["g"] = df["Geneid"].map(base)
    return dict(zip(df["g"], df[col])), {}


def read_star(path):
    """STAR ReadsPerGene.out.tab: col1 unstranded counts; first 4 rows are totals."""
    vals, special = {}, {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if f[0].startswith("N_"):
            special[f[0]] = float(f[1])
        else:
            vals[base(f[0])] = float(f[1])
    return vals, special


CONDITIONS = [
    ("STAR GeneCounts (GDC v2 rule), v36", lambda s: (f"{W}/bam/{s}_ReadsPerGene.out.tab", read_star)),
    ("htseq-count intersection-nonempty (GDC v1 rule), v36", lambda s: (f"{W}/counts/{s}_htseq_ine_v36.tsv", read_htseq)),
    ("htseq-count intersection-nonempty, v22 (v22 alignment)",
     lambda s: (f"{W}/counts/{s}_htseq_ine_v22on22.tsv", read_htseq)),
    ("htseq-count union, v36", lambda s: (f"{W}/counts/{s}_htseq_union_v36.tsv", read_htseq)),
    ("featureCounts default (unique assignment), v36", lambda s: (f"{W}/counts/{s}_fc_default_v36.tsv", read_fc)),
    ("featureCounts -O --fraction (shares ambiguous reads), v36", lambda s: (f"{W}/counts/{s}_fc_fraction_v36.tsv", read_fc)),
]

rows, extras = [], []
for s in SAMPLES:
    for label, resolve in CONDITIONS:
        path, reader = resolve(s)
        if not os.path.exists(path):
            print("missing:", path, file=sys.stderr)
            continue
        vals, special = reader(path)
        total = sum(vals.values())
        row = {"sample": s, "condition": label, "assigned_total": total}
        for gid, name in GENES.items():
            row[name] = vals.get(gid, float("nan"))
        rows.append(row)
        if special:
            extras.append({"sample": s, "condition": label,
                           **{k.lstrip("_"): v for k, v in special.items()}})

tab = pd.DataFrame(rows)
tab.to_csv(f"{W}/validation_counts.tsv", sep="\t", index=False)
pd.DataFrame(extras).to_csv(f"{W}/validation_special_counters.tsv", sep="\t", index=False)
pd.set_option("display.width", 220, "display.max_colwidth", 60)
print(tab.to_string(index=False))
print()
print(pd.DataFrame(extras).to_string(index=False))
