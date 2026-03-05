#!/usr/bin/env bash
set -euo pipefail

############################################################
# 03_pileup_phasing_numbat.sh
#
# Purpose:
#   Generate phased SNP information for downstream numbat analysis
#   using scRNA-seq BAM (Cell Ranger) and matched-normal SNP VCF.
#
# Inputs (edit paths):
#   - scRNA Cell Ranger outs:
#       <PROJECT_DIR>/scRNA/Case01/outs/
#         - possorted_genome_bam.bam
#         - filtered_feature_bc_matrix/barcodes.tsv.gz  (recommended)
#         - (or) feature_bc_matrix_cell_barcodes.csv    (if you have it)
#   - SNP VCF:
#       <PROJECT_DIR>/wgs/Case01_Normal.vcf.gz
#
# Outputs (repo-local):
#   output/numbat/phasing/Case01/
############################################################

# ==========================
# Paths (EDIT BEFORE RUN)
# ==========================
BASE_DIR=/PATH/TO/clonal_inference
DATA_DIR=${BASE_DIR}/data
OUT_DIR=${BASE_DIR}/output

# External project directory (where large inputs live)
PROJECT_DIR=/PATH/TO/PROJECT

# Reference resources (large; keep external)
GMAP=/PATH/TO/genetic_map_hg38_withX.txt.gz
PANELDIR=/PATH/TO/1000G_hg38_reference_panel

# ==========================
# Sample (EDIT)
# ==========================
SAMPLE=Case01
NCORES=24

# ==========================
# Inputs (external project)
# ==========================
CELLRANGER_OUTS=${PROJECT_DIR}/scRNA/${SAMPLE}/outs
WGS_DIR=${PROJECT_DIR}/wgs

BAM=${CELLRANGER_OUTS}/possorted_genome_bam.bam
BARCODES=${CELLRANGER_OUTS}/filtered_feature_bc_matrix/barcodes.tsv.gz
SNP_VCF=${WGS_DIR}/${SAMPLE}_Normal.vcf.gz

# ==========================
# Output (repo-local)
# ==========================
OUT_SUB=${OUT_DIR}/numbat/phasing/${SAMPLE}
mkdir -p "${OUT_SUB}"

# -----------------------------
# Find numbat script
# -----------------------------
PILEUP_AND_PHASE_R="$(Rscript -e "cat(system.file('bin/pileup_and_phase.R', package='numbat'))")"
test -f "${PILEUP_AND_PHASE_R}" || { echo "Could not locate numbat pileup_and_phase.R" >&2; exit 1; }

# ==========================
# Run
# ==========================
Rscript "${PILEUP_AND_PHASE_R}" \
  --label "${SAMPLE}" \
  --samples "${SAMPLE}" \
  --bams "${BAM}" \
  --barcodes "${BARCODES}" \
  --outdir "${OUT_SUB}" \
  --gmap "${GMAP}" \
  --snpvcf "${SNP_VCF}" \
  --paneldir "${PANELDIR}" \
  --ncores "${NCORES}"

echo "Done."
echo "  OUT: ${OUT_SUB}"






