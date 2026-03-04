#!/usr/bin/env bash
set -euo pipefail

############################################################
# numbat: pileup + phasing
#
# Purpose:
#   Generate phased SNP information for downstream numbat analysis
#   using scRNA-seq BAM (Cell Ranger) and matched-normal SNP VCF.
##
############################################################

# -----------------------------
# User-defined variables
# -----------------------------

SAMPLE="Case01"

# External project base (data lives here)
BASE_DIR="/PATH/TO/PROJECT"
CELLRANGER_OUTS="${BASE_DIR}/scRNA/${SAMPLE}/outs"
WGS_DIR="${BASE_DIR}/wgs"

# Repo root (output lives here)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
OUTDIR="${REPO_ROOT}/output/numbat/phasing/${SAMPLE}"

SNP_VCF="${WGS_DIR}/${SAMPLE}_Normal.vcf.gz"
GMAP="/PATH/TO/genetic_map_hg38_withX.txt.gz"
PANELDIR="/PATH/TO/1000G_hg38_reference_panel"

BAM="${CELLRANGER_OUTS}/cellranger_possorted_genome_bam.bam"
BARCODES="${CELLRANGER_OUTS}/feature_bc_matrix_cell_barcodes.csv"

NCORES=24

PILEUP_AND_PHASE_R="$(Rscript -e "cat(system.file('bin/pileup_and_phase.R', package='numbat'))")"

mkdir -p "${OUTDIR}"

Rscript "${PILEUP_AND_PHASE_R}" \
  --label "${SAMPLE}" \
  --samples "${SAMPLE}" \
  --bams "${BAM}" \
  --barcodes "${BARCODES}" \
  --outdir "${OUTDIR}" \
  --gmap "${GMAP}" \
  --snpvcf "${SNP_VCF}" \
  --paneldir "${PANELDIR}" \
  --ncores "${NCORES}"




