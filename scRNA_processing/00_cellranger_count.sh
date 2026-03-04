#!/usr/bin/env bash
set -euo pipefail

############################################################
# 01_cellranger_count.sh
#
# scRNA-seq preprocessing using 10x Genomics Cell Ranger
#
# Reference:
#   - refdata-gex-GRCh38-2020-A (human)
#   - refdata-gex-GRCh38-and-mm10-2020-A (PDX)
#
# Input FASTQ examples:
#
# Human tumor (Case01):
#   scrnaseq_case01_tumor_L001_R1.fastq.gz
#   scrnaseq_case01_tumor_L001_R2.fastq.gz
#   scrnaseq_case01_tumor_L002_R1.fastq.gz
#   scrnaseq_case01_tumor_L002_R2.fastq.gz
#
# PDX samples:
#   scrnaseq_pdxm_isotype_L002_R1.fastq.gz
#   scrnaseq_pdxm_isotype_L002_R2.fastq.gz
#   scrnaseq_pdxm_ly6g_L002_R1.fastq.gz
#   scrnaseq_pdxm_ly6g_L002_R2.fastq.gz
#
# Usage:
#   bash 01_cellranger_count.sh <RUN_ID> <SAMPLE_NAME> <FASTQ_DIR>
############################################################

# ==========================
# Paths (EDIT BEFORE RUN)
# ==========================
CELLRANGER=/PATH/TO/cellranger-7.0.0/cellranger

# Choose one reference:
REF=/PATH/TO/refdata-gex-GRCh38-2020-A
# REF=/PATH/TO/refdata-gex-GRCh38-and-mm10-2020-A

# ==========================
# Arguments
# ==========================
RUN_ID=$1
SAMPLE=$2
FASTQ_DIR=$3

# ==========================
# Run
# ==========================
$CELLRANGER count \
  --include-introns true \
  --id="$RUN_ID" \
  --sample="$SAMPLE" \
  --transcriptome="$REF" \
  --fastqs="$FASTQ_DIR" \
  --localcores=24 \
  --localmem=48 \
  --expect-cells=10000
