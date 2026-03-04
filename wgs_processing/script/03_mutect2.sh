#!/bin/bash
set -euo pipefail

############################################################
# 03_mutect2.sh
#
# Somatic variant calling using GATK Mutect2
# + FilterMutectCalls
#
# Input:
#   - Tumor BAM  : wgs_case01_tumor_L002_bqsr.bam
#   - Normal BAM : wgs_case01_normal_L002_bqsr.bam
#
# Output:
#   WGS_Case01.filt.vcf.gz
#
# Reference: hg38
# Sample name (SM tag): Case01
############################################################

# ==========================
# Paths (EDIT BEFORE RUN)
# ==========================
BAM_DIR=/PATH/TO/WGS_BAM_OUTPUT
OUT_DIR=/PATH/TO/MUTECT2_OUTPUT

REF_FASTA=/PATH/TO/Homo_sapiens_assembly38.fasta
GERMLINE_RESOURCE=/PATH/TO/af-only-gnomad.hg38.vcf.gz

mkdir -p "$OUT_DIR"

# ==========================
# Sample
# ==========================
SAMPLE_NAME=Case01

TUMOR_BAM=${BAM_DIR}/wgs_case01_tumor_L002_bqsr.bam
NORMAL_BAM=${BAM_DIR}/wgs_case01_normal_L002_bqsr.bam

UNFILTERED_VCF=${OUT_DIR}/${SAMPLE_NAME}.unfiltered.vcf.gz
FILTERED_VCF=${OUT_DIR}/WGS_${SAMPLE_NAME}.filt.vcf.gz

# ==========================
# 1) Mutect2
# ==========================
gatk Mutect2 \
  -R "$REF_FASTA" \
  -I "$TUMOR_BAM" \
  -I "$NORMAL_BAM" \
  -normal Case01 \
  --germline-resource "$GERMLINE_RESOURCE" \
  --native-pair-hmm-threads 24 \
  -O "$UNFILTERED_VCF"

# ==========================
# 2) FilterMutectCalls
# ==========================
gatk FilterMutectCalls \
  -R "$REF_FASTA" \
  -V "$UNFILTERED_VCF" \
  -O "$FILTERED_VCF"

# ==========================
# 3) Index
# ==========================
tabix -p vcf "$FILTERED_VCF"

echo "Mutect2 + filtering completed for ${SAMPLE_NAME}"
