#!/bin/bash
set -euo pipefail

############################################################
# 01_gatk_bwa_hg38.sh
#
# Purpose:
#   Whole-genome sequencing preprocessing (one sample per run)
#   - Alignment (BWA-MEM) to hg38
#   - Sorting and indexing (samtools)
#   - MarkDuplicates (GATK)
#   - Base Quality Score Recalibration (BQSR)
#
# Input FASTQ files (example):
#   wgs_case01_{tumor,normal}_L002_R1.fastq.gz
#   wgs_case01_{tumor,normal}_L002_R2.fastq.gz
#
############################################################

# ==========================
# Paths (EDIT BEFORE RUN)
# ==========================
BASE_DIR="/PATH/TO/wgs_processing"

DATA_DIR="${BASE_DIR}/data"
OUT_DIR="${BASE_DIR}/output"

REF_FASTA="${DATA_DIR}/Homo_sapiens_assembly38.fasta"
KNOWN_SITES_VCF="${DATA_DIR}/dbsnp138.hg38.vcf.gz"  # bgzip+tabix indexed

# ==========================
# Sample / site (EDIT)
# ==========================
SAMPLE_ID="case01"
LANE="L002"
SITE="tumor"   # set "tumor" or "normal" and run this script twice

# FASTQ prefix and read group tags
PREFIX="wgs_${SAMPLE_ID}_${SITE}_${LANE}"       # e.g., wgs_case01_tumor_L002
SM_TAG="Case01_${SITE}"                        # e.g., Case01_tumor / Case01_normal
RG_ID="${PREFIX}"

R1=${DATA_DIR}/${PREFIX}_R1.fastq.gz
R2=${DATA_DIR}/${PREFIX}_R2.fastq.gz

READGROUP="@RG\tID:${RG_ID}\tLB:lib1\tPL:ILLUMINA\tSM:${SM_TAG}\tPU:${LANE}"

# ==========================
# Output files
# ==========================
SORT_BAM="${OUT_DIR}/${PREFIX}_sort.bam"
MARKDUP_BAM="${OUT_DIR}/${PREFIX}_markdup.bam"
BQSR_BAM="${OUT_DIR}/${PREFIX}_bqsr.bam"

METRICS_TXT="${OUT_DIR}/${PREFIX}.markdup.metrics.txt"
RECAL_TXT="${OUT_DIR}/${PREFIX}_recal.txt"

# ==========================
# 1) Alignment
# ==========================
bwa mem -t 8 -K 10000000 -R "$READGROUP" "$REF_FASTA_BWA" "$R1" "$R2" \
  | samtools sort -@8 -m 10G -O bam -o "$SORT_BAM"

samtools index "$SORT_BAM"

# ==========================
# 2) MarkDuplicates
# ==========================
gatk --java-options "-Xmx50g" MarkDuplicates \
  -I "$SORT_BAM" \
  -O "$MARKDUP_BAM" \
  -M "$METRICS_TXT" \
  --TMP_DIR "$TMP_DIR"

# ==========================
# 3) BQSR
# ==========================
gatk --java-options "-Xmx50g" BaseRecalibrator \
  -I "$MARKDUP_BAM" \
  -O "$RECAL_TXT" \
  --known-sites "$KNOWN_SITES_VCF" \
  -R "$REF_FASTA_GATK"

gatk ApplyBQSR \
  -R "$REF_FASTA_GATK" \
  -I "$MARKDUP_BAM" \
  --bqsr-recal-file "$RECAL_TXT" \
  -O "$BQSR_BAM"
