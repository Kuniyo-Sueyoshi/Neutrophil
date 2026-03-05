#!/usr/bin/env python3
############################################################
# 01_convert_facets_to_phylowgs.py
#
# Convert FACETS VCF (bgzip/gzip) to PhyloWGS CNV input format.
#
# Default input (repo-local):
#   ../data/WGS.cnvfacet_Case01.vcf.gz
#
# Default output (repo-local):
#   ../output/Case01_phylowgs_cnv.tsv
#
# Optional usage:
#   python 01_convert_facets_to_phylowgs.py <input.vcf.gz> <output.tsv>
############################################################

from __future__ import annotations

import sys
import gzip
from pathlib import Path


def parse_info(info_str: str) -> dict[str, str]:
    info_dict = {}
    for item in info_str.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info_dict[key] = value
    return info_dict


def open_text_maybe_gzip(path: Path):
    # VCF from FACETS is usually bgzip/gzip; gzip module can read bgzip blocks in many cases.
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("rt")


def main():
    script_dir = Path(__file__).resolve().parent
    base_dir = script_dir.parent  # clonal_inference/

    # Defaults (repo-local)
    default_in = base_dir / "data" / "WGS.cnvfacet_Case01.vcf.gz"
    default_out = base_dir / "output" / "Case01_phylowgs_cnv.tsv"

    # Optional CLI override
    if len(sys.argv) == 1:
        input_vcf = default_in
        output_tsv = default_out
    elif len(sys.argv) == 3:
        input_vcf = Path(sys.argv[1])
        output_tsv = Path(sys.argv[2])
    else:
        print("Usage: python 01_convert_facets_to_phylowgs.py [input.vcf(.gz) output.tsv]")
        sys.exit(1)

    if not input_vcf.exists():
        raise FileNotFoundError(f"Input not found: {input_vcf}")

    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    n_written = 0
    with open_text_maybe_gzip(input_vcf) as f_in, output_tsv.open("wt") as f_out:
        f_out.write("chromosome\tstart\tend\tminor_cn\tmajor_cn\tcellular_prevalence\n")

        for line in f_in:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            chrom = fields[0].replace("chr", "")
            start = fields[1]
            info = parse_info(fields[7])

            end = info.get("END")
            tcn = info.get("TCN_EM")
            lcn = info.get("LCN_EM")
            cf  = info.get("CF_EM")

            if not (end and tcn and lcn and cf):
                continue

            try:
                tcn_i = int(round(float(tcn)))
                lcn_i = int(round(float(lcn)))
                major = tcn_i - lcn_i
                cf_f = float(cf)
            except ValueError:
                continue

            f_out.write(f"{chrom}\t{start}\t{end}\t{lcn_i}\t{major}\t{cf_f}\n")
            n_written += 1

    print("Wrote:", output_tsv)
    print("Records:", n_written)


if __name__ == "__main__":
    main()

    
