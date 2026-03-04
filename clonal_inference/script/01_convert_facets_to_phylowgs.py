#!/usr/bin/env python3

############################################################
# Convert FACETS VCF to PhyloWGS CNV input format
############################################################

import sys
import gzip

def parse_info(info_str):
    info_dict = {}
    for item in info_str.split(";"):
        if "=" in item:
            key, value = item.split("=")
            info_dict[key] = value
    return info_dict

if len(sys.argv) != 3:
    print("Usage: python convert_facets_to_phylowgs.py input.vcf.gz output.tsv")
    sys.exit(1)

# input; WGS.cnvfacet_Case01.vcf.gz
input_vcf = sys.argv[1]
output_tsv = sys.argv[2]

with gzip.open(input_vcf, 'rt') as f_in, open(output_tsv, "w") as f_out:
    f_out.write("chromosome\tstart\tend\tminor_cn\tmajor_cn\tcellular_prevalence\n")
    for line in f_in:
        if line.startswith("#"):
            continue

        fields = line.strip().split("\t")
        chrom = fields[0].replace("chr", "")
        start = fields[1]
        info = parse_info(fields[7])

        end = info.get("END")
        tcn = info.get("TCN_EM")
        lcn = info.get("LCN_EM")
        cf  = info.get("CF_EM")

        if end and tcn and lcn and cf:
            try:
                tcn = int(round(float(tcn)))
                lcn = int(round(float(lcn)))
                major = tcn - lcn
                cf_float = float(cf)

                f_out.write(f"{chrom}\t{start}\t{end}\t{lcn}\t{major}\t{cf_float}\n")

            except ValueError:
                continue
            
