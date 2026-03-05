#!/usr/bin/env Rscript

############################################################
# 02_tximport_summarize.R
#
# Purpose:
#   Summarize Salmon transcript quantification (quant.sf)
#   to gene-level expression using tximport.
#
# Inputs (repo-local):
#   output/<DATASET>/salmon/<sample>/quant.sf
#   data/tx2gene.tsv
#
# Output (repo-local):
#   output/<DATASET>/gene_expression.rds
############################################################

suppressPackageStartupMessages({
  library(tximport)
  library(readr)
  library(dplyr)
})

# ==========================
# Paths (EDIT BEFORE RUN)
# ==========================
BASE_DIR <- "/PATH/TO/bulkRNA_processing"

DATA_DIR <- file.path(BASE_DIR, "data")
OUT_DIR  <- file.path(BASE_DIR, "output")

# Dataset (EDIT)
DATASET <- "CellM_Histamine"   # or "Neutrophil_PDX"

SALMON_DIR <- file.path(OUT_DIR, DATASET, "salmon")
TX2GENE_TSV <- file.path(DATA_DIR, "tx2gene.tsv")

OUT_RDS <- file.path(OUT_DIR, DATASET, "gene_expression.rds")

# -----------------------------
# quant.sf files
# -----------------------------
files <- list.files(
  SALMON_DIR,
  recursive = TRUE,
  pattern = "quant\\.sf$",
  full.names = TRUE
)

sample_names <- basename(dirname(files))
names(files) <- sample_names

print(files)

# -----------------------------
# tximport (transcript-level)
# -----------------------------
tx.exp <- tximport(
  files,
  type = "salmon",
  txOut = TRUE
)

# -----------------------------
# Tx2gene
# -----------------------------
tx2gene <- read_tsv(
  TX2GENE_TSV,
  show_col_types = FALSE
)

# -----------------------------
# Gene-level summarization
# -----------------------------
gene.exp <- summarizeToGene(
  tx.exp,
  tx2gene,
  countsFromAbundance = "lengthScaledTPM"
)

# -----------------------------
# Save
# -----------------------------
dir.create(dirname(OUT_RDS), showWarnings = FALSE, recursive = TRUE)

saveRDS(gene.exp, OUT_RDS)

