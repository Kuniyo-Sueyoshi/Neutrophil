#!/usr/bin/env Rscript
############################################################
# 00_tximport_summarize.R
#
# Summarize Salmon quantification (quant.sf) to gene-level
# using tximport.
#
# Usage:
#   Rscript 00_tximport_summarize.R <DATASET> <TX2GENE_TSV>
#
# DATASET:
#   CellM_Histamine
#   Neutrophil_PDX
#
# Expected input structure:
#   sf/<DATASET>/<SAMPLE>/quant.sf
#
# Output:
#   output/bulkRNAseq_<DATASET>.rds
############################################################

suppressPackageStartupMessages({
  library(tximport)
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

DATASET     <- args[1]
TX2GENE_TSV <- args[2]

# -----------------------------
# Directory settings
# -----------------------------
BASE_DIR <- ".."  # bulkRNA_processing/

SF_DIR  <- file.path(BASE_DIR, "sf", DATASET)
OUT_DIR <- file.path(BASE_DIR, "output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OUT_RDS <- file.path(OUT_DIR, paste0("bulkRNAseq_", DATASET, ".rds"))

# -----------------------------
# quant.sf files
# -----------------------------
files <- list.files(SF_DIR, recursive = TRUE,
                    pattern = "quant\\.sf$",
                    full.names = TRUE)

sample_names <- basename(dirname(files))
names(files) <- sample_names

print(files)

# -----------------------------
# tximport (transcript-level)
# -----------------------------
tx.exp <- tximport(files, type = "salmon", txOut = TRUE)

# -----------------------------
# Tx2gene (assumed 2-column table)
# Columns:
#   TXNAME
#   GENE
# -----------------------------
tx2gene <- read_tsv(TX2GENE_TSV, show_col_types = FALSE)

# -----------------------------
# Gene-level summarization
# -----------------------------
gene.exp <- summarizeToGene(
  tx.exp,
  tx2gene,
  countsFromAbundance = "lengthScaledTPM"
)

saveRDS(gene.exp, OUT_RDS)
