#!/usr/bin/env bash
set -euo pipefail

############################################################
# 01_salmon_quant.sh
#
# Bulk RNA-seq transcript quantification using Salmon
#
# Usage:
#   bash 01_salmon_quant.sh <DATASET>
#
# DATASET options:
#   CellM_Histamine
#   Neutrophil_PDX
#
# FASTQ naming convention (paired-end, trimmed):
#   *_1_trim.fq.gz
#   *_2_trim.fq.gz
#
# Example output directories:
#
# CellM_Histamine:
#   DW_1/   DW_2/   DW_3/
#   HA_1/   HA_2/   HA_3/
#
# Neutrophil_PDX:
#   Neu_N24/ Neu_N25/ Neu_N26/
#   Neu_S13/ Neu_S14/ Neu_S15/
############################################################

# ==========================
# Paths (EDIT BEFORE RUN)
# ==========================
BASE_DIR=/PATH/TO/BULK_PROJECT

DATASET=$1

FASTQ_DIR=${BASE_DIR}/fastq/${DATASET}
OUT_DIR=${BASE_DIR}/sf/${DATASET}

# Choose reference index:
REF_INDEX=/PATH/TO/SalmonhgTx
# REF_INDEX=/PATH/TO/SalmonmmTx

THREADS=16

mkdir -p "${OUT_DIR}"

# ==========================
# Run Salmon
# ==========================
for fq1 in ${FASTQ_DIR}/*_1_trim.fq.gz
do
    fq2=${fq1/_1_trim.fq.gz/_2_trim.fq.gz}
    sample=$(basename "${fq1}" _1_trim.fq.gz)
    sample_out=${OUT_DIR}/${sample}

    salmon quant \
        -i "${REF_INDEX}" \
        -l A \
        -1 "${fq1}" \
        -2 "${fq2}" \
        -p ${THREADS} \
        -o "${sample_out}" \
        --validateMappings
done
