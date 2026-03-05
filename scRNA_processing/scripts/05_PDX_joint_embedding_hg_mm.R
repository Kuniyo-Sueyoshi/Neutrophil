#!/usr/bin/env Rscript
############################################################
# 05_PDX_joint_embedding_hg_mm.R
#
# Purpose
#   Joint embedding of:
#     - human tumor cells (hg)
#     - mouse stromal cells (mm → converted to human symbols)
#
# Steps
#   1) Load annotated hg tumor object
#   2) Load annotated mm stroma object
#   3) Convert mouse gene symbols → human symbols
#   4) Collapse duplicated genes (per-cell MAX)
#   5) Merge hg tumor + converted mm stroma
#   6) SCTransform + Harmony integration
#   7) Save integrated Seurat object
#
# Inputs (repo-local)
#   output/PDX_integration_mm/PDX_mm_Integrated.rds
#   output/PDX_integration_hg/PDX_hg_Integrated.rds
#   data/HMD_HumanPhenotype.rpt
#
# Output
#   output/PDX_joint_embedding_hg_mm/scRNAseq_PDX_NeuDep.rds
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(readr)
})

set.seed(123)

# ==========================
# Paths (EDIT BEFORE RUN)
# ==========================
BASE_DIR <- "/PATH/TO/scRNA_processing"

DATA_DIR <- file.path(BASE_DIR, "data")
OUT_BASE <- file.path(BASE_DIR, "output")
R_DIR    <- file.path(BASE_DIR, "R")

OUT_DIR <- file.path(OUT_BASE, "PDX_joint_embedding_hg_mm")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ==========================
# Inputs
# ==========================
mouse_rds <- file.path(
  OUT_BASE,
  "PDX_integration_mm",
  "PDX_mm_Integrated.rds"
)

tumor_rds <- file.path(
  OUT_BASE,
  "PDX_integration_hg",
  "PDX_hg_Integrated.rds"
)

homolog_tsv <- file.path(
  DATA_DIR,
  "HMD_HumanPhenotype.rpt"
)

out_rds <- file.path(
  OUT_DIR,
  "scRNAseq_PDX_NeuDep.rds"
)

# ==========================
# Helper functions
# ==========================
source(file.path(R_DIR, "fn.cluster.SCT.R"))

# -----------------------------
# Sanity checks
# -----------------------------
stopifnot(file.exists(mouse_rds))
stopifnot(file.exists(tumor_rds))
stopifnot(file.exists(homolog_tsv))

# ============================================================
# 1) Load objects
# ============================================================
d.mouse <- readRDS(mouse_rds)
d.meso  <- readRDS(tumor_rds)

# Ensure tumor label
d.meso$cell_middle <- "Tumor"

# ============================================================
# 2) Mouse → human gene symbol conversion
# ============================================================
homolog <- read_tsv(
  homolog_tsv,
  col_names = c(
    "Symbol.hg",
    "key",
    "Symbol.mm",
    "MGI.ID",
    "Phenotype",
    "X6"
  ),
  show_col_types = FALSE
)

mm2hg <- homolog$Symbol.hg
names(mm2hg) <- homolog$Symbol.mm

mm_counts <- GetAssayData(
  d.mouse,
  assay = "RNA",
  slot = "counts"
)

mm_genes <- rownames(mm_counts)
hg_genes <- unname(mm2hg[mm_genes])

keep <- !is.na(hg_genes) & hg_genes != "NULL"

mm_counts <- mm_counts[keep, , drop = FALSE]
hg_genes  <- hg_genes[keep]

# ============================================================
# 3) Collapse duplicated genes (sparse-safe MAX)
# ============================================================
collapse_by_max_sparse <- function(mat, new_rownames) {
  stopifnot(inherits(mat, "dgCMatrix"))
  stopifnot(nrow(mat) == length(new_rownames))

  idx_list <- split(seq_along(new_rownames), new_rownames)
  uniq <- names(idx_list)

  out <- Matrix(0, nrow = length(uniq), ncol = ncol(mat), sparse = TRUE)
  rownames(out) <- uniq
  colnames(out) <- colnames(mat)

  for (i in seq_along(uniq)) {
    ridx <- idx_list[[i]]
    if (length(ridx) == 1L) {
      out[i, ] <- mat[ridx, ]
    } else {
      v <- mat[ridx[1], ]
      for (k in ridx[-1]) v <- pmax(v, mat[k, ])
      out[i, ] <- v
    }
  }
  out
}

hg_counts <- collapse_by_max_sparse(mm_counts, hg_genes)

# ============================================================
# 4) Create converted mouse object
# ============================================================
d.mouse.hg <- CreateSeuratObject(
  counts    = hg_counts,
  meta.data = d.mouse@meta.data,
  project   = "mouse_as_human_symbols"
)

# ============================================================
# 5) Merge tumor + stroma
# ============================================================
data <- merge(d.meso, d.mouse.hg)

# ============================================================
# 6) SCT + Harmony integration
# ============================================================
data <- fn.cluster.SCT(
  data = data,
  batch = "batch",
  vars.to.regress = c(
    "G2M.Score",
    "S.Score",
    "percent.mt"
  )
)

# ============================================================
# 7) Coarse label
# ============================================================
if (!"cell_coarse" %in% colnames(data@meta.data)) {
  data$cell_coarse <- NA_character_
}
data$cell_coarse[data$cell_middle == "Tumor"] <- "Tumor"

# ============================================================
# 8) Save
# ============================================================
saveRDS(data, out_rds)
