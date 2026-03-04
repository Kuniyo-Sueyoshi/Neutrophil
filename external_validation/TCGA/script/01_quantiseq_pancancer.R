#!/usr/bin/env Rscript

############################################################
# 01_quantiseq_pancancer.R
#
# Purpose:
#   Run quanTIseq (TIL10) deconvolution on TCGA PanCan RNAseqV2
#   expression matrix (logRSEM-like values) after per-sample
#   per-million scaling.
#
# Inputs (generalized):
#   - data/TCGA_RNAseqV2_logRSEM.tsv
#     * rows: gene symbols (column "gene_id")
#     * cols: TCGA sample barcodes (e.g., TCGA-XX-YYYY-...)
#
# Outputs:
#   - output/Immune10_deconvolv_TCGA.PanCan.tsv
#   - output/Immune10_deconvolv_TCGA.PanCan.annotated.tsv
#
# Notes:
#   - The input file is described by GDC PanCanAtlas as batch-corrected,
#     length-corrected, isoform-aware RSEM outputs. This script treats it
#     as an expression matrix and rescales each sample to 1e6 total.
#   - Sample barcodes are shortened to patient-level (TCGA-XX-YYYY).
############################################################


suppressPackageStartupMessages({
  library(tidyverse)
  library(quantiseqr)
})

set.seed(123)

# -----------------------------
# Paths (EDIT BEFORE RUN)
# -----------------------------
BASE_DIR <- "/PATH/TO/external_validation/TCGA"

IN_TSV  <- file.path(BASE_DIR, "data", "TCGA_RNAseqV2_logRSEM.tsv")
OUT_DIR <- file.path(BASE_DIR, "output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
IN_CLIN_TSV <- file.path(BASE_DIR, "data", "TCGA_clinical_annotation_PanCan.tsv")

OUT_TSV <- file.path(OUT_DIR, "Immune10_deconvolv_TCGA.PanCan.tsv")
OUT_TSV_ANNOT <- file.path(OUT_DIR, "Immune10_deconvolv_TCGA.PanCan.annotated.tsv")

# -----------------------------
# 1) Load expression matrix
# -----------------------------
df_RSEM <- read_tsv(IN_TSV, show_col_types = FALSE) %>%
  column_to_rownames("gene_id")

# -----------------------------
# 2) Standardize sample IDs to patient-level (TCGA-XX-YYYY)
# -----------------------------
colnames(df_RSEM) <- colnames(df_RSEM) %>%
  str_split(pattern = "-", simplify = FALSE) %>%
  sapply(function(x) paste(x[1], x[2], x[3], sep = "-"))

# -----------------------------
# 3) Per-million scaling (TPM-like)
# -----------------------------
df_TPM <- apply(df_RSEM, 2, function(x) {
  x <- as.numeric(x)
  x / sum(x) * 1e6
})

# keep rownames
rownames(df_TPM) <- rownames(df_RSEM)

# -----------------------------
# 4) quanTIseq (TIL10)
# -----------------------------
ti_tcga <- quantiseqr::run_quantiseq(
  expression_data = df_TPM,
  signature_matrix = "TIL10",
  is_arraydata = FALSE,
  is_tumordata = TRUE,
  scale_mRNA = TRUE
)

write_tsv(ti_tcga, OUT_TSV)

# -----------------------------
# 5) Merge with clinical annotations (optional)
# -----------------------------
meta.data <- read_tsv(IN_CLIN_TSV, show_col_types = FALSE) %>%
  # Expect at least: tumor_barcode, bcr_patient_barcode, acronym,
  # days_to_last_followup, days_to_death, vital_status
  mutate(
    days_to_death = na_if(days_to_death, "[Not Applicable]"),
    days_to_last_followup = na_if(days_to_last_followup, "[Not Applicable]"),
    days_to_censor = ifelse(!is.na(days_to_death), days_to_death, days_to_last_followup),
    days_to_censor = as.numeric(days_to_censor),
    status = ifelse(vital_status == "Dead", 1, 0)
  ) %>%
  select(
    bcr_patient_barcode,
    acronym,
    days_to_last_followup,
    days_to_death,
    vital_status,
    status,
    days_to_censor
  )

ti_tcga_annotated <- merge(
  ti_tcga,
  meta.data,
  by.x = "Sample",
  by.y = "bcr_patient_barcode",
  all.x = TRUE
)

write_tsv(ti_tcga_annotated, OUT_TSV_ANNOT)
