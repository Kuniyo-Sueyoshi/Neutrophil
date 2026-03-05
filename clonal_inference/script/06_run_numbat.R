#!/usr/bin/env Rscript

############################################################
# 06_run_numbat.R
#
# Purpose:
#   Infer tumor subclones from scRNA-seq data using numbat.
#
# Inputs:
#   1) Integrated Seurat object (hg38)                 [external]
#   2) Allele counts from pileup_and_phase.R           [repo-local output]
#   3) WGS-derived coarse CNV segment file             [repo-local output]
#   4) Reference cell-type expression matrix           [external or repo data]
#
# Outputs (repo-local):
#   output/numbat/results/<SAMPLE>/
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyverse)
  library(numbat)
})

set.seed(123)

# ==========================
# Paths (EDIT BEFORE RUN)
# ==========================
BASE_DIR <- "/PATH/TO/clonal_inference"
DATA_DIR <- file.path(BASE_DIR, "data")
OUT_DIR  <- file.path(BASE_DIR, "output")

# Large external inputs live here
PROJECT_DIR <- "/PATH/TO/PROJECT"

# ==========================
# Sample
# ==========================
SAMPLE <- "Case01"

# ==========================
# Inputs
# ==========================
SEURAT_RDS <- file.path(PROJECT_DIR, "data", paste0(SAMPLE, "_seurat_integrated.rds"))

ALLELE_TSV <- file.path(
  OUT_DIR, "numbat", "phasing", SAMPLE,
  paste0(SAMPLE, "_allele_counts.tsv.gz")
)

SEG_FILE <- file.path(
  OUT_DIR, "numbat", "segments", SAMPLE,
  paste0(SAMPLE, "_segments_simplified.tsv")
)

# Reference profile (external recommended if large)
REF_CELLTYPE_RDS <- file.path(
  PROJECT_DIR, "numbat", "reference",
  "ref_celltype_expression.rds"
)

# ==========================
# Output
# ==========================
OUT_NUMBAT <- file.path(OUT_DIR, "numbat", "results", SAMPLE)
dir.create(OUT_NUMBAT, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) Load Seurat object
# -----------------------------
data <- readRDS(SEURAT_RDS)
data_tumor <- subset(data, subset = cell_coarse == "Mesothelial")

count_mat <- GetAssayData(data_tumor, slot = "counts", assay = "RNA")

# -----------------------------
# 2) Load reference cell-type profile
# -----------------------------
ref_celltype <- readRDS(REF_CELLTYPE_RDS)

# -----------------------------
# 3) Load allele counts
# -----------------------------
df_allele <- read.delim(ALLELE_TSV, sep = "\t", stringsAsFactors = FALSE)

# -----------------------------
# 4) Load CNV segments
# -----------------------------
df_segs_consensus <- read.delim(SEG_FILE, sep = "\t", stringsAsFactors = FALSE)

df_segs_loh <- df_segs_consensus %>%
  filter(cnv_state == "loh") %>%
  mutate(CHROM = as.factor(CHROM), loh = TRUE)

# -----------------------------
# 5) Detect clonal LOH (bulk-level)
# -----------------------------
psdbulk <- numbat::get_bulk(
  count_mat   = count_mat,
  lambdas_ref = ref_celltype,
  df_allele   = df_allele,
  gtf         = numbat::gtf_hg38,
  segs_loh    = NULL  # or df_segs_loh
)

clonal_loh <- numbat::detect_clonal_loh(bulk = psdbulk)

# -----------------------------
# 6) Run numbat
# -----------------------------
out <- run_numbat(
  count_mat = count_mat,
  ref_celltype = ref_celltype,
  df_allele = df_allele,
  genome = "hg38",
  t = 1e-5,				     # Defalt
  ncores = 16,
  # segs_consensus_fix = df_segs_consensus,  # optional
  call_clonal_loh = TRUE,
  plot = TRUE,
  out_dir = OUT_DIR
)

saveRDS(out, file.path(OUT_NUMBAT, "numbat_results.rds"))



