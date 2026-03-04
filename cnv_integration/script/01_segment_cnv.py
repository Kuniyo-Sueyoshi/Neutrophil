#!/usr/bin/env python3

############################################################
# 01_segment_cnv.py
#
# Convert FACETS VCF to window-segmented CNV table
# (per-sample; indexed by Segmented_Chr)
#
# Input:
#   - FACETS VCF(s): set via FACETS_DIR and FILE_GLOB
#   - Segmented_Chr template: ../output/chr_segmented.tsv
#
# Output:
#   ../output/facetTCN_SegmentedChr_<SAMPLE>.tsv
############################################################

import gzip
import re
import glob
from pathlib import Path

import numpy as np
import pandas as pd

from utils_make_SegChr import make_SegChr


def parse_info_field(info_str: str) -> dict:
    """Parse VCF INFO field into dict."""
    out = {}
    for item in info_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
    return out


def main():
    script_dir = Path(__file__).resolve().parent
    base_dir = script_dir.parent  # cnv_integration/

    # -----------------------------
    # Paths (EDIT BEFORE RUN)
    # -----------------------------
    FACETS_DIR = Path("/PATH/TO/cnv_facets")  # e.g., /data/share/WGS/.../cnv_facets/
    FILE_GLOB = "M03_2_T.vcf.gz"              # e.g., "*.vcf.gz" or "M03_2_T.vcf.gz"

    out_dir = base_dir / "output"
    out_dir.mkdir(parents=True, exist_ok=True)

    tmp_dir = out_dir / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    template_file = out_dir / "chr_segmented.tsv"  # produced by 00_segment_template.py

    in_files = sorted(glob.glob(str(FACETS_DIR / FILE_GLOB)))
    print(in_files)

    for in_file in in_files:
        in_path = Path(in_file)

        # sample name from filename: "<sample>.vcf.gz"
        sample = re.sub(r"\.vcf\.gz$", "", in_path.name)

        out_file = out_dir / f"facetTCN_SegmentedChr_{sample}.tsv"
        tmp_file = tmp_dir / f"tmp_facetTCN_SegmentedChr_{sample}.tsv"

        print("Sample:", sample)

        ploidy = None
        n = 0

        with gzip.open(in_path, "rt") as fr:
            for line in fr:
                if line.startswith("##ploidy="):
                    ploidy = round(float(line.strip().replace("##ploidy=", "")))
                    continue

                if line.startswith("#"):
                    continue

                items = line.rstrip("\n").split("\t")
                chrom = items[0]
                pos = int(items[1])
                flt = items[6]
                info = parse_info_field(items[7])

                # required keys
                if "END" not in info or "TCN_EM" not in info:
                    continue

                if ploidy is None:
                    # FACETS VCF should contain ploidy header; if not, skip sample
                    continue

                tcn = int(round(float(info["TCN_EM"])))
                end = int(info["END"])

                # original logic
                mTCN = int(tcn - ploidy) if tcn == ploidy else int((tcn - ploidy) / abs(tcn - ploidy))

                # filters (keep original intent)
                if flt != "PASS":
                    continue
                if tcn == ploidy:
                    continue
                if chrom == "chrX":
                    continue

                chr_num = int(chrom.replace("chr", ""))

                chr_segs, chr_seg_lens = make_SegChr(CHR=chr_num, start=pos, end=end)

                d_seg_tcn = pd.DataFrame(
                    {
                        "Segmented_Chr": np.array(chr_segs, dtype=int),
                        "TCN": np.array([tcn] * len(chr_segs)),
                        "mTCN": np.array([mTCN] * len(chr_segs)),
                        "Segmented_Chr_Len": np.array(chr_seg_lens, dtype=float),
                    }
                )

                if n == 0:
                    d_seg_tcn.to_csv(tmp_file, sep="\t", index=False)
                else:
                    d_seg_tcn.to_csv(tmp_file, sep="\t", mode="a", index=False, header=False)
                n += 1

        # Summarize duplicated Segmented_Chr
        d = pd.read_csv(tmp_file, sep="\t")
        d["TCN"] = d["TCN"] * d["Segmented_Chr_Len"]
        d["mTCN"] = d["mTCN"] * d["Segmented_Chr_Len"]

        dg = d.groupby("Segmented_Chr")[["TCN", "mTCN"]].sum().round({"TCN": 3, "mTCN": 3})
        print(dg.head())

        # merge with empty template
        d_seg_all = pd.read_csv(template_file, sep="\t").set_index("Segmented_Chr")
        dg_all = pd.merge(d_seg_all, dg, how="outer", left_index=True, right_index=True, sort=False)

        dg_all.fillna({"TCN": ploidy, "mTCN": 0}).to_csv(out_file, sep="\t", index=True)
        print("Wrote:", out_file)


if __name__ == "__main__":
    main()

