#!/usr/bin/env bash
set -euo pipefail

############################################################
# 02_prepare_and_run_phylowgs.sh
#
# Purpose:
#   Filter Mutect2 VCF for PhyloWGS input and run PhyloWGS.
#
# Inputs (repo-local):
#   data/WGS_Case01.filt.vcf.gz               # Mutect2 filtered VCF (bgzip+tabix recommended)
#   output/Case01_phylowgs_cnv.tsv            # CNV TSV (from 01_convert_facets_to_phylowgs.py)
#
# Outputs (repo-local):
#   output/phylowgs/Case01/
#     - Case01.filtered.phylowgs.vcf
#     - input/ (ssm_data.txt, cnv_data.txt, ...)
#     - chains/ (trees.zip, ...)
#     - Case01.summary.json.gz, Case01.mutations.json.gz, Case01.mutass.zip
############################################################

# ==========================
# Paths (EDIT BEFORE RUN)
# ==========================
BASE_DIR=/PATH/TO/clonal_inference

DATA_DIR=${BASE_DIR}/data
OUT_BASE=${BASE_DIR}/output

# Tools (edit to your environment)
PHYLOWGS_DIR=/PATH/TO/phylowgs
PHYLOWGS_PARSER=${PHYLOWGS_DIR}/parser/create_phylowgs_inputs.py
PHYLOWGS_RUN=${PHYLOWGS_DIR}/multievolve.py
WRITE_RESULTS=${PHYLOWGS_DIR}/write_results.py

# Dependencies assumed on PATH:
#   bcftools
#   python2  (PhyloWGS often requires python2; keep as-is if your install does)
PYTHON2=python2

# ==========================
# Sample (EDIT)
# ==========================
CASE_ID=Case01

# ==========================
# Inputs
# ==========================
INPUT_VCF=${DATA_DIR}/WGS_${CASE_ID}.filt.vcf.gz
CNV_TSV=${OUT_BASE}/${CASE_ID}_phylowgs_cnv.tsv

# ==========================
# Outputs
# ==========================
OUT_DIR=${OUT_BASE}/phylowgs/${CASE_ID}
mkdir -p "${OUT_DIR}"

FILTERED_VCF=${OUT_DIR}/${CASE_ID}.filtered.phylowgs.vcf

# PhyloWGS output locations (explicit; no cd)
INPUT_DIR=${OUT_DIR}/input
CHAINS_DIR=${OUT_DIR}/chains

SUMMARY_JSON=${OUT_DIR}/${CASE_ID}.summary.json.gz
MUT_JSON=${OUT_DIR}/${CASE_ID}.mutations.json.gz
MUTASS_ZIP=${OUT_DIR}/${CASE_ID}.mutass.zip

# ==========================================================
# 1) Filter VCF for PhyloWGS
# ==========================================================
bcftools view \
  -i 'INFO/TLOD > 10 && INFO/NLOD > 3 && INFO/DP >= 50' \
  "${INPUT_VCF}" \
  > "${FILTERED_VCF}"

# ==========================================================
# 2) Create PhyloWGS inputs (writes to OUT_DIR/input)
# ==========================================================
mkdir -p "${INPUT_DIR}"

${PYTHON2} "${PHYLOWGS_PARSER}" \
  --cnvs "${CASE_ID}=${CNV_TSV}" \
  --vcf-type "${CASE_ID}=mutect_smchet" \
  "${CASE_ID}=${FILTERED_VCF}" \
  --tumor-sample "${CASE_ID}" \
  --output-dir "${INPUT_DIR}"

# ==========================================================
# 3) Run PhyloWGS (writes to OUT_DIR/chains)
# ==========================================================
mkdir -p "${CHAINS_DIR}"

${PYTHON2} "${PHYLOWGS_RUN}" \
  --num-chains 16 \
  --ssms "${INPUT_DIR}/ssm_data.txt" \
  --cnvs "${INPUT_DIR}/cnv_data.txt" \
  --burnin-samples 500 \
  --mcmc-samples 2000 \
  --output-dir "${CHAINS_DIR}"

# ==========================================================
# 4) Summarize results
# ==========================================================
${PYTHON2} "${WRITE_RESULTS}" \
  "${CASE_ID}" \
  "${CHAINS_DIR}/trees.zip" \
  "${SUMMARY_JSON}" \
  "${MUT_JSON}" \
  "${MUTASS_ZIP}"

echo "Done."
echo "  OUT_DIR: ${OUT_DIR}"


