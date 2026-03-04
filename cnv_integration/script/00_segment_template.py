#!/usr/bin/env python3

############################################################
# 00_segment_template.py
#
# Create a Segmented_Chr index table from chromosome lengths.
#
# Input:
#   ../data/chr_len.tsv
#
# Output:
#   ../output/chr_segmented.tsv
############################################################

from pathlib import Path
import pandas as pd

# If utils_make_SegChr.py is in the same directory as this script:
from utils_make_SegChr import make_SegChr


def main():
    script_dir = Path(__file__).resolve().parent
    base_dir = script_dir.parent  # cnv_integration/

    in_file = base_dir / "data" / "chr_len.tsv"
    out_dir = base_dir / "output"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "chr_segmented.tsv"

    chrseg_ls = []

    with in_file.open("rt") as fr:
        next(fr)  # skip header
        for line in fr:
            items = line.rstrip("\n").split("\t")
            chr_id = int(items[0])
            chr_len = int(items[1])

            chr_segs, _chr_seg_lens = make_SegChr(CHR=chr_id, start=1, end=chr_len)
            chrseg_ls.extend(chr_segs)

    d_seg_all = pd.DataFrame(index=chrseg_ls)
    d_seg_all.index.name = "Segmented_Chr"
    d_seg_all.to_csv(out_file, sep="\t")

    print("Wrote:", out_file)


if __name__ == "__main__":
    main()

