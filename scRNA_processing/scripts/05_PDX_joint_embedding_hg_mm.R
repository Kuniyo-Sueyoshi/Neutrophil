#!/usr/bin/env Rscript
# ============================================================
# 05_PDX_joint_embedding_hg_mm.R
# Joint embedding (hg tumor + mm stroma)
#   - Load hg tumor Seurat object (annotated)
#   - Load mm stroma Seurat object (annotated)
#   - Convert mouse genes -> human symbols (HomoloGene table)
#     * drop unmapped genes
#     * collapse duplicated human symbols by per-cell MAX (sparse-safe)
#   - Merge hg tumor + mouse(stroma, human-symbolized)
#   - SCTransform (glmGamPoi) + Harmony integration (fn.cluster.SCT)
#   - Save integrated Seurat object
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(readr)
})

set.seed(123)

# -----------------------------
# 0) Resolve repo root (stable) + helpers
# -----------------------------
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
SCRIPT_DIR <- normalizePath(dirname(script_path))
REPO_ROOT  <- normalizePath(file.path(SCRIPT_DIR, ".."))  # scRNA_processing/

source(file.path(REPO_ROOT, "R", "fn.cluster.SCT.R"))  # fn.cluster.SCT()

# -----------------------------
# 1) Inputs (generalized)
# -----------------------------
BASE_DIR <- "/PATH/TO/PROJECT"          # EDIT
IN_DIR   <- file.path(BASE_DIR, "data") # where *.rds and reference live

# Output (repo-local)
OUT_DIR <- file.path(REPO_ROOT, "output", "PDX_joint_embedding_hg_mm")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

mouse_rds   <- file.path(IN_DIR, "PDX_mm_integrated.rds")   # mm stroma (annotated)
tumor_rds   <- file.path(IN_DIR, "PDX_hg_Integrated.rds")   # hg tumor (annotated)
homolog_tsv <- file.path(IN_DIR, "HMD_HumanPhenotype.rpt")  # EDIT location/name

out_rds <- file.path(OUT_DIR, "scRNAseq_PDX_NeuDep.rds")

# -----------------------------
# 2) Load objects
# -----------------------------
d.mouse <- readRDS(mouse_rds)
d.meso  <- readRDS(tumor_rds)

# ensure tumor label exists (methods explicitly states tumor vs stroma)
d.meso$cell_middle <- "Tumor"

# -----------------------------
# 3) Mouse -> human gene symbol conversion (HomoloGene)
# -----------------------------
# accept possible 6th column (often present as X6)
homolog <- read_tsv(
  homolog_tsv,
  col_names = c("Symbol.hg", "key", "Symbol.mm", "MGI.ID", "Phenotype", "X6"),
  show_col_types = FALSE
)

# named character vector: mouse_symbol -> human_symbol
mm2hg <- homolog$Symbol.hg
names(mm2hg) <- homolog$Symbol.mm

mm_counts <- GetAssayData(d.mouse, assay = "RNA", slot = "counts")
mm_genes  <- rownames(mm_counts)

hg_genes <- unname(mm2hg[mm_genes])

keep <- !is.na(hg_genes) & hg_genes != "NULL"
mm_counts <- mm_counts[keep, , drop = FALSE]
hg_genes  <- hg_genes[keep]

# ---- collapse duplicated human symbols by MAX per cell (sparse-safe) ----
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

d.mouse.hg <- CreateSeuratObject(
  counts    = hg_counts,
  meta.data = d.mouse@meta.data,
  project   = "mouse_as_human_symbols"
)

# -----------------------------
# 4) Merge + SCT/Harmony joint embedding
# -----------------------------
data <- merge(d.meso, d.mouse.hg)

data <- fn.cluster.SCT(
  data = data,
  batch = "batch",
  vars.to.regress = c("G2M.Score", "S.Score", "percent.mt")
)

# -----------------------------
# 5) Coarse label: Tumor vs non-tumor
# -----------------------------
if (!"cell_coarse" %in% colnames(data@meta.data)) {
  data$cell_coarse <- NA_character_
}
data$cell_coarse[data$cell_middle == "Tumor"] <- "Tumor"

# -----------------------------
# 6) Save
# -----------------------------
saveRDS(data, out_rds)
