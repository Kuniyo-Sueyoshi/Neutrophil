#!/usr/bin/env Rscript

# ============================================================
# PDX scRNA-seq preprocessing (per-sample)
# ============================================================
# Sample:
sampleID <- "PDX.M_Ly6G"
# sampleID <- "PDX.M_Isotype"

# Overview
# - Load filtered 10x matrix
# - Ambient RNA correction evaluated (SoupX); not applied if low-complexity
# - Standard QC filtering
# - Doublet removal (scDblFinder)
# - Cross-species multiplet exclusion using hg/(hg+mm) fraction
# - Split into human (hg) and mouse (mm) compartments
# - Save Seurat objects for downstream integration
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(scDblFinder)
  library(SingleCellExperiment)
})

# -----------------------------
# Sample
# -----------------------------
sampleID <- "PDX.M_Ly6G"
# sampleID <- "PDX.M_Isotype"

# -----------------------------
# Paths (generalized)
# -----------------------------
BASE_DIR <- "/PATH/TO/10X/"

# repo-local output (scRNA_processing/output/<sampleID>)
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
SCRIPT_DIR <- normalizePath(dirname(script_path))
REPO_ROOT  <- normalizePath(file.path(SCRIPT_DIR, ".."))  # scRNA_processing
OUT_DIR    <- file.path(REPO_ROOT, "output", "PDX_split", sampleID)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# QC thresholds
# -----------------------------
UMI.cutoff     <- 500
Feature.cutoff <- 200
Mt.cutoff      <- 15
Hb.cutoff      <- 5

hg_dom  <- 0.8
mm_dom  <- 0.2
dbl_margin <- 0.2  # (kept as-is; unused)

# -----------------------------
# Helper functions
# -----------------------------
source(file.path(REPO_ROOT, "R", "fn.QCmetrics.R"))

# -----------------------------
# 1) Load 10X matrix
# -----------------------------
counts <- Read10X(
  data.dir = file.path(BASE_DIR, sampleID, "outs", "filtered_feature_bc_matrix")
)
obj <- CreateSeuratObject(counts = counts, min.cells = 1, min.features = 1)

# -----------------------------
# 3) QC metrics & filtering
# -----------------------------
obj <- add_qc_metrics(obj, cap = TRUE)

obj <- subset(
  obj,
  subset =
    nCount_RNA > UMI.cutoff &
    nFeature_RNA > Feature.cutoff &
    percent.mt < Mt.cutoff &
    percent.hb < Hb.cutoff
)

# -----------------------------
# 4) Doublet detection
# -----------------------------
sce <- as.SingleCellExperiment(obj)
sce <- scDblFinder(sce)
obj$scDblFinder.class <- colData(sce)$scDblFinder.class

# -----------------------------
# 5) Cross-species read fraction
# -----------------------------
gene.species <- rownames(obj) |>
  gsub("GRCh38.*", "GRCh38", x = _) |>
  gsub("mm10.*", "mm10", x = _)

counts.mat <- GetAssayData(obj, slot = "counts")

species.sum <- apply(counts.mat, 2, function(vec)
  tapply(vec, INDEX = gene.species, sum)
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

# -----------------------------
# 6) Split hg / mm
# -----------------------------
cell.species <- ifelse(
  obj$ratio.hg_mm >= hg_dom, "hg",
  ifelse(obj$ratio.hg_mm <= mm_dom, "mm", NA)
)

obj$species <- cell.species

obj.hg <- subset(obj, subset = species == "hg")
obj.mm <- subset(obj, subset = species == "mm")

RenameGenes <- function(x, pattern) {
  newnames <- gsub(pattern, "", rownames(x))
  rownames(x) <- newnames
  x
}

obj.hg <- RenameGenes(obj.hg, "GRCh38-")
obj.mm <- RenameGenes(obj.mm, "mm10-")

# -----------------------------
# 7) Save
# -----------------------------
saveRDS(obj.hg, file = file.path(OUT_DIR, paste0(sampleID, "_hg.rds")))
saveRDS(obj.mm, file = file.path(OUT_DIR, paste0(sampleID, "_mm.rds")))
