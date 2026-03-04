#!/usr/bin/env python3

############################################################
# Merge adjacent CNV segments with identical cnv_state
############################################################

#!/usr/bin/env python3
import argparse
import pandas as pd
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--seg", required=True)
parser.add_argument("--sample", required=True)
parser.add_argument("--outdir", required=True)
args = parser.parse_args()

SEG_FILE = args.seg
SAMPLE = args.sample
OUTDIR = args.outdir

os.makedirs(OUTDIR, exist_ok=True)

df = pd.read_csv(SEG_FILE, sep="\t")
if df.empty:
    out_file = os.path.join(OUTDIR, f"{SAMPLE}_segments_simplified.tsv")
    # write header-only file to avoid crash and make downstream behavior explicit
    pd.DataFrame(columns=["CHROM", "seg", "seg_start", "seg_end", "cnv_state"]).to_csv(
        out_file, sep="\t", index=False
    )
    print("Simplified segment file written to:", out_file, file=sys.stderr)
    sys.exit(0)

df = df.sort_values(by=["CHROM", "seg_start"]).reset_index(drop=True)

merged = []
current = df.iloc[0].to_dict()

for i in range(1, len(df)):
    row = df.iloc[i]
    if (row["CHROM"] == current["CHROM"] and row["cnv_state"] == current["cnv_state"]):
        current["seg_end"] = max(current["seg_end"], row["seg_end"])
    else:
        merged.append(current)
        current = row.to_dict()

merged.append(current)
result = pd.DataFrame(merged)

# reassign segment labels (BUGFIX: remove 26-limit)
result["seg"] = result.groupby("CHROM").cumcount().apply(lambda x: f"seg{x:04d}")
result["seg"] = result["CHROM"].astype(str) + "_" + result["seg"]

out_file = os.path.join(OUTDIR, f"{SAMPLE}_segments_simplified.tsv")
result[["CHROM", "seg", "seg_start", "seg_end", "cnv_state"]].to_csv(
    out_file, sep="\t", index=False
)

print("Simplified segment file written to:", out_file)

