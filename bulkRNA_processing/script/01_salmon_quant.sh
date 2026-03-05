#!/usr/bin/env bash
set -euo pipefail

############################################################
# 00_salmon_quant.sh
#
# Purpose:
#   Bulk RNA-seq transcript quantification using Salmon.
#
# Inputs (repo-local):
#   data/<DATASET>/*_1_trim.fq.gz
#   data/<DATASET>/*_2_trim.fq.gz
#
# Output (repo-local):
#   output/<DATASET>/salmon/<sample>/quant.sf
#
# Notes:
#   - Paired-end reads, trimmed FASTQs with suffix *_1_trim.fq.gz / *_2_trim.fq.gz
############################################################

# ==========================
# Paths (EDIT BEFORE RUN)
# ==========================
BASE_DIR="/PATH/TO/bulkRNA_processing"

DATA_DIR="${BASE_DIR}/data"
OUT_DIR="${BASE_DIR}/output"

# Dataset selector (EDIT)
DATASET="CellM_Histamine"   # or "Neutrophil_PDX"

FASTQ_DIR="${DATA_DIR}/${DATASET}"
SALMON_OUT="${OUT_DIR}/${DATASET}/salmon"

# Salmon index (EDIT)
REF_INDEX="/PATH/TO/SalmonIndex"   # e.g., hgTx or mmTx

THREADS=16
mkdir -p "$SALMON_OUT"

# ==========================
# Run Salmon
# ==========================

for fq1 in "${FASTQ_DIR}"/*_1_trim.fq.gz; do
  fq2="${fq1/_1_trim.fq.gz/_2_trim.fq.gz}"
  sample="$(basename "$fq1" _1_trim.fq.gz)"
  sample_out="${SALMON_OUT}/${sample}"

  salmon quant \
    -i "$REF_INDEX" \
    -l A \
    -1 "$fq1" \
    -2 "$fq2" \
    -p "$THREADS" \
    -o "$sample_out" \
    --validateMappings
done

