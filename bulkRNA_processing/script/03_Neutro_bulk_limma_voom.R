#!/usr/bin/env Rscript
############################################################
### About input file: bulkRNAseq_Neutrophil.rds
# Bulk RNA-seq of mouse neutrophils sorted from PDX tumors.
# Quantification:
#   Salmon quantification followed by xenograft-aware reference assignment
#   (human vs mouse; mouse-mapped transcripts retained), then summarized to
#   gene-level with tximport.
#
# Input object (tximport gene-level, saved as RDS) contains at least:
#   $counts    : raw gene counts (genes x samples) for limma-voom
#   $abundance : abundance values from tximport (genes x samples; used for GSVA)

### About this script
# Differential expression and pathway scoring:
#   - limma-voom for DEG analysis
#   - GSVA using curated neutrophil function gene sets

# Contrast:
#   M vs K (M - K)
#   M: PDX.M-derived neutrophils (N24–N26)
#   K: PDX.K-derived neutrophils (S13–S15)

# Outputs:
#   1) DEG table (TSV)
#   2) GSVA scores (TSV; all gene sets)
#   3) GSVA selected pathways (TSV; subset)
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
IN_RDS   <- file.path(REPO_ROOT, "data", "bulkRNAseq_Neutrophil.rds")
GSET_TSV <- file.path(REPO_ROOT, "data", "neutrophil_functions.txt")

OUT_DIR <- file.path(REPO_ROOT, "output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OUT_DEG_TSV          <- file.path(OUT_DIR, "DEG_M_vs_K.tsv")
OUT_GSVA_ALL_TSV     <- file.path(OUT_DIR, "GSVA_all_pathways.tsv")
OUT_GSVA_SELECT_TSV  <- file.path(OUT_DIR, "GSVA_selected_pathways.tsv")

# -----------------------------
# 1) Load tximport output
# -----------------------------
data <- readRDS(IN_RDS)

counts <- data$counts
abund  <- data$abundance

stopifnot(is.matrix(counts) || inherits(counts, "Matrix"))
stopifnot(is.matrix(abund)  || inherits(abund,  "Matrix"))
stopifnot(all(colnames(counts) == colnames(abund)))

# -----------------------------
# 2) Design matrix (M vs K)
# -----------------------------
sampleTable <- data.frame(
  condition = factor(c(rep("M", 3), rep("K", 3)))
)
rownames(sampleTable) <- colnames(counts)

design <- model.matrix(~ 0 + condition, data = sampleTable)
colnames(design) <- gsub("condition", "", colnames(design))

# -----------------------------
# 3) edgeR object + simple filtering
# -----------------------------
y <- DGEList(counts)

cutoff <- 10
drop <- which(apply(y$counts, 1, max) < cutoff)
if (length(drop) > 0) {
  y <- y[-drop, , keep.lib.sizes = FALSE]
}

y <- calcNormFactors(y)

# -----------------------------
# 4) limma-voom + DE
# -----------------------------
v <- voom(y, design)

fit <- lmFit(v, design)
cont.mat <- makeContrasts(MvsK = M - K, levels = design)
fit <- contrasts.fit(fit, cont.mat)
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
# -----------------------------
df_gset <- read_tsv(GSET_TSV, show_col_types = FALSE)

if (!all(c("Pathway", "Genes") %in% colnames(df_gset))) {
  colnames(df_gset)[1:2] <- c("Pathway", "Genes")
}

gset_list <- setNames(
  lapply(df_gset$Genes, function(x) str_split(x, pattern = ",\\s*")[[1]]),
  df_gset$Pathway
)

# -----------------------------
# 6) GSVA on log(TPM+1)
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

write_tsv(gsva_df, OUT_GSVA_ALL_TSV)

# -----------------------------
# 7) Export selected pathways
# -----------------------------
selected <- c("Angiogenesis", "ECM remodeling", "Interferon signaling")

gsva_selected <- gsva_df %>%
  filter(Pathway %in% selected)

col_rename <- c(
  "T_Neu_N24" = "TAN_M.rep1",
  "T_Neu_N25" = "TAN_M.rep2",
  "T_Neu_N26" = "TAN_M.rep3",
  "T_Neu_S13" = "TAN_K.rep1",
  "T_Neu_S14" = "TAN_K.rep2",
  "T_Neu_S15" = "TAN_K.rep3"
)
common_cols <- intersect(names(col_rename), colnames(gsva_selected))
if (length(common_cols) > 0) {
  gsva_selected <- gsva_selected %>%
    rename(!!!col_rename[common_cols])
}

write_tsv(gsva_selected, OUT_GSVA_SELECT_TSV)
