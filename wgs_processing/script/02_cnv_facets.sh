#!/bin/bash
set -euo pipefail

############################################################
# 02_cnv_facets.sh
#
# Allele-specific copy number analysis using FACETS
#
# Input:
#   - Tumor BAM  : wgs_case01_tumor_L002_bqsr.bam
#   - Normal BAM : wgs_case01_normal_L002_bqsr.bam
#
# Output:
#   FACETS results (segmentation, purity/ploidy estimates)
#
# Sample ID: Case01
############################################################

# ==========================
# Paths (EDIT BEFORE RUN)
# ==========================
BAM_DIR="/PATH/TO/WGS_BAM_OUTPUT"

REF_FASTA="/PATH/TO/Homo_sapiens_assembly38.fasta"
KNOWN_SITES="/PATH/TO/dbsnp_138.hg38.vcf.gz"
GENE_BED="/PATH/TO/hg38_gene.bed"

# ==========================
# Sample
# ==========================
SAMPLE_NAME="Case01"

TUMOR_BAM="${BAM_DIR}/wgs_case01_tumor_L002_bqsr.bam"
NORMAL_BAM="${BAM_DIR}/wgs_case01_normal_L002_bqsr.bam"

# ==========================
# Output directories (repo-local)
# ==========================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

OUT_DIR="${REPO_ROOT}/output/${SAMPLE_NAME}/facets"
mkdir -p "$OUT_DIR"

# ==========================
# Run FACETS
# ==========================
cnv_facets.R \
  -t "$TUMOR_BAM" \
  -n "$NORMAL_BAM" \
  -vcf "$KNOWN_SITES" \
  --annotation "$GENE_BED" \
  -o "${OUT_DIR}/${SAMPLE_NAME}"
