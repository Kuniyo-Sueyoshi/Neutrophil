#!/usr/bin/env python3

############################################################
# Convert FACETS VCF to coarse CNV segment file (for numbat)
#
# Input:
#   FACETS VCF (tumor sample; hg38)
#
# Output:
#   seg file:
#   CHROM | seg | seg_start | seg_end | cnv_state
#
# Usage:
#   python wgs_facets_to_seg.py \
#       --vcf input.vcf.gz \
#       --sample Case01 \
#       --outdir ./output
############################################################

import argparse
import pandas as pd
import gzip
import re
from collections import defaultdict
import string
import os

# -----------------------------
# Argument parser
# -----------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--vcf", required=True)
parser.add_argument("--sample", required=True)
parser.add_argument("--outdir", required=True)
args = parser.parse_args()

VCF_FILE = args.vcf
SAMPLE = args.sample
OUTDIR = args.outdir

os.makedirs(OUTDIR, exist_ok=True)

# -----------------------------
# SVTYPE → cnv_state mapping
# -----------------------------
svtype_to_cnv_state = {
    "DEL": "del",
    "DUP": "amp",
    "DUP-LOH": "amp",
    "HEMIZYG": "loh",
    "LOH": "loh",
    "NEUTR": "neu"
}

# -----------------------------
# Read VCF
# -----------------------------
records = []

open_func = gzip.open if VCF_FILE.endswith(".gz") else open

with open_func(VCF_FILE, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue

        parts = line.strip().split("\t")
        chrom_raw, pos, _, _, _, _, _, info = parts[:8]

        # chr1 → 1
        chrom_match = re.match(r"chr(\d+|X|Y)", chrom_raw)
        if not chrom_match:
            continue

        chrom = chrom_match.group(1)

        # convert X/Y to numeric labels if needed
        if chrom == "X":
            chrom = 23
        elif chrom == "Y":
            chrom = 24
        else:
            chrom = int(chrom)

        start = int(pos)

        info_dict = dict(
            item.split("=")
            for item in info.split(";")
            if "=" in item
        )

        end = int(info_dict.get("END", start))
        svtype = info_dict.get("SVTYPE", "NEUTR")

        cnv_state = svtype_to_cnv_state.get(svtype, "neu")

        records.append({
            "CHROM": chrom,
            "seg_start": start,
            "seg_end": end,
            "cnv_state": cnv_state
        })

df = pd.DataFrame(records)

# -----------------------------
# Assign segment IDs
# -----------------------------
letters = list(string.ascii_lowercase)
double_letters = [a + b for a in letters for b in letters]

seg_ids = []
chrom_seg_index = defaultdict(int)

for _, row in df.iterrows():
    chrom = row["CHROM"]
    idx = chrom_seg_index[chrom]
    seg_label = f"{chrom}{double_letters[idx]}"
    seg_ids.append(seg_label)
    chrom_seg_index[chrom] += 1

df["seg"] = seg_ids

df_final = df[["CHROM", "seg", "seg_start", "seg_end", "cnv_state"]]

out_file = os.path.join(OUTDIR, f"{SAMPLE}_segments.tsv")
df_final.to_csv(out_file, sep="\t", index=False)

print("Segment file written to:", out_file)
