#!/usr/bin/env Rscript

############################################################
# 03_Neutro_bulk_limma_voom.R
#
# Purpose:
#   Bulk RNA-seq analysis of mouse neutrophils sorted from PDX tumors:
#     - limma-voom for DEG (M vs K)
#     - GSVA using curated neutrophil function gene sets
#
# Input (repo-local):
#   output/Neutrophil_PDX/gene_expression.rds
#     - tximport output containing at least:
#         $counts    (genes x samples)
#         $abundance (genes x samples)
#   data/neutrophil_functions.tsv
#     - two columns: Pathway, Genes (comma-separated)
#
# Outputs (repo-local):
#   output/Neutrophil_PDX/DEG_M_vs_K.tsv
#   output/Neutrophil_PDX/GSVA_pathways.tsv
#
# Notes:
#   - Assumes 6 samples in the order: M (n=3), K (n=3)
############################################################

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(tibble)
  library(dplyr)
  library(readr)
  library(stringr)
  library(GSVA)
})

# ==========================
# Paths (EDIT BEFORE RUN)
# ==========================
BASE_DIR <- "/PATH/TO/bulkRNA_processing"

DATA_DIR <- file.path(BASE_DIR, "data")
OUT_BASE <- file.path(BASE_DIR, "output")

DATASET <- "Neutrophil_PDX"

IN_RDS   <- file.path(OUT_BASE, DATASET, "gene_expression.rds")
GSET_TSV <- file.path(DATA_DIR, "neutrophil_functions.tsv")  # Pathway, Genes

OUT_DIR <- file.path(OUT_BASE, DATASET)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OUT_DEG_TSV     <- file.path(OUT_DIR, "DEG_M_vs_K.tsv")
OUT_GSVA_TSV    <- file.path(OUT_DIR, "GSVA_pathways.tsv")

# -----------------------------
# 1) Load tximport output
# -----------------------------
stopifnot(file.exists(IN_RDS))
obj <- readRDS(IN_RDS)

stopifnot(!is.null(obj$counts), !is.null(obj$abundance))
counts <- obj$counts
abund  <- obj$abundance

stopifnot(all(colnames(counts) == colnames(abund)))

# -----------------------------
# 2) Design matrix (M vs K)
# -----------------------------
stopifnot(ncol(counts) == 6)

sample_table <- data.frame(
  condition = factor(c(rep("M", 3), rep("K", 3)), levels = c("K", "M")),
  row.names = colnames(counts)
)

design <- model.matrix(~ 0 + condition, data = sample_table)
colnames(design) <- gsub("^condition", "", colnames(design))

# -----------------------------
# 3) edgeR: filtering + TMM
# -----------------------------
y <- DGEList(counts)
keep <- filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

# -----------------------------
# 4) limma-voom + DE (M - K)
# -----------------------------
v <- voom(y, design)

fit <- lmFit(v, design)
cont <- makeContrasts(MvsK = M - K, levels = design)
fit <- contrasts.fit(fit, cont)
fit <- eBayes(fit)

deg <- topTable(
  fit,
  coef = "MvsK",
  adjust = "fdr",
  number = Inf,
  sort.by = "P"
) %>%
  rownames_to_column("Gene") %>%
  filter(!grepl("description", Gene))

write_tsv(deg, OUT_DEG_TSV)

# -----------------------------
# 5) Gene sets (neutrophil functions)
# - "Angiogenesis"
# - "ECM remodeling"
# - "Interferon signaling"
# (See supplementary Tables)
# -----------------------------
df_gset <- read_tsv(GSET_TSV, show_col_types = FALSE)

# allow flexible header naming (first two cols)
if (!all(c("Pathway", "Genes") %in% colnames(df_gset))) {
  colnames(df_gset)[1:2] <- c("Pathway", "Genes")
}

gset_list <- setNames(
  lapply(df_gset$Genes, function(x) str_split(x, pattern = ",\\s*")[[1]]),
  df_gset$Pathway
)

# -----------------------------
# 6) GSVA on log(TPM+1) computed from abundance
# -----------------------------
keep_abund <- apply(abund, 1, max) != 0
tpm <- abund[keep_abund, , drop = FALSE]
tpm <- apply(tpm, 2, function(x) x * 1e6 / sum(x))

gsva_mat <- GSVA::gsva(
  expr = log(tpm + 1),
  gset.idx.list = gset_list,
  kcdf = "Gaussian"
)

gsva_df <- as.data.frame(gsva_mat) %>%
  rownames_to_column("Pathway")

write_tsv(gsva_df, OUT_GSVA_TSV)




