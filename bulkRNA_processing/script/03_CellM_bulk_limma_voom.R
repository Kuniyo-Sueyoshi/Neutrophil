#!/usr/bin/env Rscript

############################################################
# 03_CellM_bulk_limma_voom.R
#
# Purpose:
#   Bulk RNA-seq differential expression analysis (CellM):
#     - edgeR (filtering + TMM normalization)
#     - limma-voom linear modeling
#
# Input (repo-local):
#   output/CellM_Histamine/gene_expression.rds   (tximport output; counts slot required)
#
# Output (repo-local):
#   output/CellM_Histamine/DEG_CellM_HA_vs_DW.tsv
#
# Notes:
#   - Assumes 6 samples in the order: DW (n=3), HA (n=3)
############################################################

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(tibble)
  library(dplyr)
})

set.seed(123)

# ==========================
# Paths (EDIT BEFORE RUN)
# ==========================
BASE_DIR <- "/PATH/TO/bulkRNA_processing"

DATASET <- "CellM_Histamine"

OUT_DIR <- file.path(BASE_DIR, "output", DATASET)
IN_RDS  <- file.path(OUT_DIR, "gene_expression.rds")

OUT_TSV <- file.path(OUT_DIR, "DEG_CellM_HA_vs_DW.tsv")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1) Load gene-level counts
# -----------------------------
stopifnot(file.exists(IN_RDS))
obj <- readRDS(IN_RDS)

stopifnot(!is.null(obj$counts))
counts <- obj$counts

# -----------------------------
# 2) Experimental design
# -----------------------------
stopifnot(ncol(counts) == 6)

sample_table <- data.frame(
  condition = factor(c(rep("DW", 3), rep("HA", 3)), levels = c("DW", "HA")),
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
# 4) limma-voom + contrast
# -----------------------------
v <- voom(y, design)
fit <- lmFit(v, design)

cont <- makeContrasts(HAvsDW = HA - DW, levels = design)
fit <- contrasts.fit(fit, cont)
fit <- eBayes(fit)

# -----------------------------
# 5) Export DEG table
# -----------------------------
deg <- topTable(
  fit,
  coef = "HAvsDW",
  adjust = "fdr",
  number = Inf,
  sort.by = "P"
) %>%
  rownames_to_column("Gene") %>%
  filter(!grepl("description", Gene))

write.table(
  deg,
  OUT_TSV,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
