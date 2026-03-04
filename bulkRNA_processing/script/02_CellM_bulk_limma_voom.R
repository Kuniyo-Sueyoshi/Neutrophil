#!/usr/bin/env Rscript
############################################################
# Bulk RNA-seq differential expression analysis
# Sample: Cell.M
#
# Pipeline:
#   Salmon quantification
#   → tximport gene-level summarization
#   → edgeR (TMM normalization)
#   → limma-voom linear modeling
#
# Input: bulkRNAseq_cellM_Histamine.rds
#   gene-level count matrix (tximport output; counts slot required)
#
# Output:
#   DEG table (HA vs DW)
############################################################

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(tibble)
  library(dplyr)
})

set.seed(123)

# -----------------------------
# Resolve repo root (stable)
# -----------------------------
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
SCRIPT_DIR <- normalizePath(dirname(script_path))
REPO_ROOT  <- normalizePath(file.path(SCRIPT_DIR, ".."))  # bulkRNA_processing/

# -----------------------------
# Paths (repo-local)
# -----------------------------
IN_RDS  <- file.path(REPO_ROOT, "data", "bulkRNAseq_cellM_Histamine.rds")
OUT_DIR <- file.path(REPO_ROOT, "output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OUT_TSV <- file.path(OUT_DIR, "DEG_CellM_HA_vs_DW.tsv")

# -----------------------------
# 1) Load gene-level counts
# -----------------------------
data <- readRDS(IN_RDS)
counts <- data$counts

# -----------------------------
# 2) Experimental design
# -----------------------------
sampleTable <- data.frame(
  condition = factor(c(rep("DW", 3), rep("HA", 3)))
)
rownames(sampleTable) <- colnames(counts)

design <- model.matrix(~ 0 + condition, data = sampleTable)
colnames(design) <- gsub("condition", "", colnames(design))

# -----------------------------
# 3) edgeR object + filtering
# -----------------------------
y <- DGEList(counts)
keep <- filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

# -----------------------------
# 4) limma-voom
# -----------------------------
v <- voom(y, design)
fit <- lmFit(v, design)

cont.mat <- makeContrasts(HAvsDW = HA - DW, levels = design)
fit <- contrasts.fit(fit, cont.mat)
fit <- eBayes(fit)

# -----------------------------
# 5) Extract DEG table
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
