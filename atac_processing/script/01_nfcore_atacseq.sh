#!/usr/bin/env bash
set -euo pipefail

############################################################
# 01_nfcore_atacseq.sh
#
# Purpose:
#   Process ATAC-seq data using nf-core/atacseq (Nextflow).
#
# Inputs (repo-local):
#   data/samplesheet.csv
#
# Output (repo-local):
#   output/   (nf-core/atacseq results)
#
# Notes:
#   - Genome: GRCh38
#   - Read length: 150 bp
############################################################

# ==========================
# Paths (EDIT BEFORE RUN)
# ==========================
BASE_DIR="/PATH/TO/atac_processing"

DATA_DIR="${BASE_DIR}/data"
OUT_DIR="${BASE_DIR}/output"

SAMPLESHEET="${DATA_DIR}/samplesheet.csv"

# Nextflow temp directory (optional)
NXF_TMPDIR="/PATH/TO/nextflow_tmp"

# nf-core options
GENOME="GRCh38"
READ_LENGTH=150
PROFILE="docker"   # or "singularity", "conda" depending on your environment

mkdir -p "$OUT_DIR"

# ==========================
# Run nf-core/atacseq
# ==========================
export NXF_OPTS="-Djava.io.tmpdir=${NXF_TMPDIR}"

nextflow run nf-core/atacseq \
  --input "$SAMPLESHEET" \
  --outdir "$OUT_DIR" \
  --genome "$GENOME" \
  --read_length "$READ_LENGTH" \
  -profile "$PROFILE"
