#!/usr/bin/env Rscript

############################################################
# 02_PDX_split_hg_mm.R
#
# Purpose:
#   Preprocess PDX scRNA-seq (per-sample; Isotype or Ly6G):
#     - Load filtered 10x matrix
#     - QC filtering
#     - Doublet removal (scDblFinder)
#     - Cross-species multiplet exclusion using hg/(hg+mm)
#     - Split into human (hg) and mouse (mm)
#
# Inputs:
#   <TENX_DIR>/<sampleID>/outs/filtered_feature_bc_matrix/
#
# Outputs (repo-local):
#   output/PDX_split/<sampleID>/<sampleID>_hg.rds
#   output/PDX_split/<sampleID>/<sampleID>_mm.rds
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(scDblFinder)
  library(SingleCellExperiment)
})

set.seed(123)

# ==========================
# Paths (EDIT BEFORE RUN)
# ==========================
BASE_DIR <- "/PATH/TO/scRNA_processing"

DATA_DIR <- file.path(BASE_DIR, "data")
OUT_BASE <- file.path(BASE_DIR, "output")
R_DIR    <- file.path(BASE_DIR, "R")

# External 10x directory (EDIT)
TENX_DIR <- "/PATH/TO/10X"

# ==========================
# Sample (EDIT)
# ==========================
sampleID <- "PDX.M_Ly6G"
# sampleID <- "PDX.M_Isotype"

OUT_DIR <- file.path(OUT_BASE, "PDX_split", sampleID)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ==========================
# QC thresholds
# ==========================
UMI.cutoff     <- 500
Feature.cutoff <- 200
Mt.cutoff      <- 15
Hb.cutoff      <- 5

hg_dom <- 0.8
mm_dom <- 0.2

# ==========================
# Helper functions
# ==========================
source(file.path(R_DIR, "fn.QCmetrics.R"))

RenameGenes <- function(x, pattern) {
  newnames <- gsub(pattern, "", rownames(x))
  rownames(x) <- newnames
  x
}

# ============================================================
# 1) Load 10X matrix
# ============================================================
counts <- Read10X(
  data.dir = file.path(TENX_DIR, sampleID, "outs", "filtered_feature_bc_matrix")
)

obj <- CreateSeuratObject(counts = counts, min.cells = 1, min.features = 1)
obj$batch <- sampleID

# ============================================================
# 2) QC metrics & filtering
# ============================================================
obj <- add_qc_metrics(obj, cap = TRUE)

obj <- subset(
  obj,
  subset =
    nCount_RNA   > UMI.cutoff &
    nFeature_RNA > Feature.cutoff &
    percent.mt   < Mt.cutoff &
    percent.hb   < Hb.cutoff
)

# ============================================================
# 3) Doublet detection
# ============================================================
sce <- as.SingleCellExperiment(obj)
sce <- scDblFinder(sce)

obj$scDblFinder.class <- colData(sce)$scDblFinder.class

# ============================================================
# 4) Cross-species read fraction
# ============================================================
gene.species <- rownames(obj) |>
  gsub("GRCh38.*", "GRCh38", x = _) |>
  gsub("mm10.*", "mm10", x = _)

counts.mat <- GetAssayData(obj, slot = "counts")

species.sum <- apply(
  counts.mat,
  2,
  function(vec) tapply(vec, INDEX = gene.species, sum)
)

ratio.hg_mm <- apply(species.sum, 2, function(vec) {
  total <- sum(vec)
  hg <- if ("GRCh38" %in% names(vec)) vec[["GRCh38"]] else 0
  if (total == 0) return(NA_real_)
  hg / total
})

obj$ratio.hg_mm <- ratio.hg_mm

keep <- (
  obj$scDblFinder.class == "singlet" &
  !(ratio.hg_mm > mm_dom & ratio.hg_mm < hg_dom)
)

obj <- obj[, keep]

# ============================================================
# 5) Split hg / mm
# ============================================================
obj$species <- ifelse(
  obj$ratio.hg_mm >= hg_dom, "hg",
  ifelse(obj$ratio.hg_mm <= mm_dom, "mm", NA)
)

obj.hg <- subset(obj, subset = species == "hg")
obj.mm <- subset(obj, subset = species == "mm")

obj.hg <- RenameGenes(obj.hg, "GRCh38-")
obj.mm <- RenameGenes(obj.mm, "mm10-")

# ============================================================
# 6) Save
# ============================================================
saveRDS(obj.hg, file = file.path(OUT_DIR, paste0(sampleID, "_hg.rds")))
saveRDS(obj.mm, file = file.path(OUT_DIR, paste0(sampleID, "_mm.rds")))

