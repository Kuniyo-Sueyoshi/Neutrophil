#!/bin/bash
set -euo pipefail

############################################################
# 01_gatk_bwa_hg38.sh
#
# Whole-genome sequencing preprocessing pipeline
# - Alignment (BWA-MEM) to hg38
# - Sorting and indexing (samtools)
# - MarkDuplicates (GATK)
# - Base Quality Score Recalibration (BQSR)
#
# Input FASTQ files:
#   wgs_case01_normal_L002_R1.fastq.gz
#   wgs_case01_normal_L002_R2.fastq.gz
#   wgs_case01_tumor_L002_R1.fastq.gz
#   wgs_case01_tumor_L002_R2.fastq.gz
#
# Sample name (SM tag): Case01
############################################################

# ==========================
# Paths (EDIT BEFORE RUN)
# ==========================
BAM_DIR="/PATH/TO/WGS_BAM_OUTPUT"

REF_FASTA="/PATH/TO/Homo_sapiens_assembly38.fasta"
GERMLINE_RESOURCE="/PATH/TO/af-only-gnomad.hg38.vcf.gz"

# ==========================
# Sample
# ==========================
SAMPLE_NAME="Case01"

TUMOR_BAM="${BAM_DIR}/wgs_case01_tumor_L002_bqsr.bam"
NORMAL_BAM="${BAM_DIR}/wgs_case01_normal_L002_bqsr.bam"

# IMPORTANT: Mutect2 -normal expects the NORMAL sample name (SM tag) in the BAM header.
# If your normal SM tag is the same as SAMPLE_NAME, keep as-is.
NORMAL_SM="${SAMPLE_NAME}"

# ==========================
# Output directories (repo-local)
# ==========================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

OUT_DIR="${REPO_ROOT}/output/${SAMPLE_NAME}/mutect2"
mkdir -p "$OUT_DIR"

UNFILTERED_VCF="${OUT_DIR}/${SAMPLE_NAME}.unfiltered.vcf.gz"
FILTERED_VCF="${OUT_DIR}/WGS_${SAMPLE_NAME}.filt.vcf.gz"

# ==========================
# 1) Mutect2
# ==========================
gatk Mutect2 \
  -R "$REF_FASTA" \
  -I "$TUMOR_BAM" \
  -I "$NORMAL_BAM" \
  -normal "$NORMAL_SM" \
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

