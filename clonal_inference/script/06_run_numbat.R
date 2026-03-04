#!/usr/bin/env Rscript

############################################################
# numbat clonal inference 
#
# Purpose:
#   Infer tumor subclones from scRNA-seq data using numbat.
#
# Inputs:
#   1) Integrated Seurat object (hg38)
#   2) Allele counts from pileup_and_phase.R
#   3) WGS-derived coarse CNV segment file
#   4) Reference cell-type expression matrix
#
# Notes:
#   Reference cell-type profiles were generated from
#   independent scRNA-seq datasets (Case12, 13, 15)
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyverse)
  library(numbat)
})

set.seed(123)

SAMPLE <- "Case01"

# repo root = one level above this script directory
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
SCRIPT_DIR <- normalizePath(dirname(script_path))
REPO_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."))

# External project base (large inputs live here)
BASE_DIR <- "/PATH/TO/PROJECT"

SEURAT_RDS <- file.path(BASE_DIR, "data", paste0(SAMPLE, "_seurat_integrated.rds"))

# Outputs from this repo's pipeline (now standardized)
ALLELE_TSV <- file.path(REPO_ROOT, "output", "numbat", "phasing", SAMPLE,
                        paste0(SAMPLE, "_allele_counts.tsv.gz"))

SEG_FILE <- file.path(REPO_ROOT, "output", "numbat", "seg", SAMPLE,
                      paste0(SAMPLE, "_segments_simplified.tsv"))

# Reference profile (choose: keep in repo data/ or external BASE_DIR)
REF_CELLTYPE_RDS <- file.path(BASE_DIR, "clonal_inference", "numbat", "reference",
                              "ref_celltype_expression.rds")

OUT_DIR <- file.path(REPO_ROOT, "output", "numbat", "results", SAMPLE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# 1) Load Seurat object
data <- readRDS(SEURAT_RDS)

data_tumor <- subset(data, subset = cell_coarse == "Mesothelial")

count_mat <- GetAssayData(
  data_tumor,
  slot = "counts",
  assay = "RNA"
)

# 2) Load reference cell-type profile
ref_celltype <- readRDS(REF_CELLTYPE_RDS)

# 3) Load allele counts
df_allele <- read.delim(ALLELE_TSV, sep = "\t", stringsAsFactors = FALSE)

# 4) WGS-derived CNV segments
df_segs_consensus <- read.delim(SEG_FILE, sep = "\t", stringsAsFactors = FALSE)

df_segs_loh <- df_segs_consensus %>%
  filter(cnv_state == "loh") %>%
  mutate(CHROM = as.factor(CHROM)) %>%
  mutate(loh = TRUE)

# 5) Detect clonal LOH (bulk-level)
psdbulk <- numbat::get_bulk(
  count_mat = count_mat,
  lambdas_ref = ref_celltype,
  df_allele = df_allele,
  gtf = numbat::gtf_hg38,
  segs_loh = NULL   # or df_segs_loh
)

clonal_loh <- numbat::detect_clonal_loh(bulk = psdbulk)

# 6) Run numbat
out <- run_numbat(
  count_mat = count_mat,
  ref_celltype = ref_celltype,
  df_allele = df_allele,
  genome = "hg38",
  t = 1e-5,
  ncores = 16,
  # segs_consensus_fix = df_segs_consensus,  # optional
  call_clonal_loh = TRUE,
  plot = TRUE,
  out_dir = OUT_DIR
)

saveRDS(out, file.path(OUT_DIR, "numbat_results.rds"))

