#!/usr/bin/env Rscript
# ============================================================
# 04_PDX_sample_integration_mm.R
# PDX (mouse; mm) sample integration: PDX.M.Iso vs PDX.M.Ly6G
#
# Steps
#   1) Load + merge mm Seurat objects
#   2) QC metrics annotation
#   3) Cell-cycle scoring (human Seurat cc.genes -> mouse orthologs via gorth)
#   4) SCTransform (glmGamPoi) + Harmony integration + clustering
#   5) Remove ambient-like / low-quality clusters (explicit IDs) + re-run
#   6) Subclustering for ambiguous populations (myeloid boundary; lymphoid contamination)
#   7) Final cell-type labels (fine / middle / coarse)
#   8) Save integrated object
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(gprofiler2)
})

set.seed(123)

# -----------------------------
# 0) Parameters (generalized)
# -----------------------------
sampleIDs <- c("PDX.M.Iso", "PDX.M.Ly6G")

# Resolve repo root (scRNA_processing/)
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
SCRIPT_DIR <- normalizePath(dirname(script_path))
REPO_ROOT  <- normalizePath(file.path(SCRIPT_DIR, ".."))

# Inputs (keep external if you want; here keep as-is but make paths explicit)
BASE_DIR <- "/PATH/TO/PROJECT"
IN_DIR   <- file.path(BASE_DIR, "data")

# Output (repo-local)
OUT_DIR  <- file.path(REPO_ROOT, "output", "PDX_integration_mm")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Input naming convention: {sampleID}_mm.rds
in_files <- file.path(IN_DIR, paste0(sampleIDs, "_mm.rds"))

omit_clusters_res02 <- c("5", "13")
refine_clusters_1 <- c("1", "4")
refine_clusters_2 <- c("3", "12")

# -----------------------------
# Helper functions (repo-local)
# -----------------------------
source(file.path(REPO_ROOT, "R", "fn.cluster.SCT.R"))    # fn.cluster.SCT()
source(file.path(REPO_ROOT, "R", "fn.QCmetrics.R"))      # add_qc_metrics()

# -----------------------------
# 1) Load and merge mm objects
# -----------------------------
data <- merge(readRDS(in_files[1]), readRDS(in_files[2]))

# Require batch in meta.data (as you intended)
if (!"batch" %in% colnames(data@meta.data)) {
  stop("meta.data$batch is missing in the input objects. Please ensure *_mm.rds include a 'batch' column.")
}

data$condition <- ifelse(grepl("Isotype", data$batch, ignore.case = TRUE), "Isotype", "Ly6G")

# -----------------------------
# 2) QC metrics annotation
# -----------------------------
data <- add_qc_metrics(data, cap = TRUE)

# -----------------------------
# 3) Cell-cycle scoring (mouse orthologs)
# -----------------------------
data("cc.genes", package = "Seurat")

mm.s.genes <- gprofiler2::gorth(
  query = cc.genes$s.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name

mm.g2m.genes <- gprofiler2::gorth(
  query = cc.genes$g2m.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name

data <- Seurat::CellCycleScoring(
  object = data,
  s.features = mm.s.genes,
  g2m.features = mm.g2m.genes,
  verbose = FALSE
)

# -----------------------------
# 4) Integration (SCT + Harmony) and clustering
# -----------------------------
data <- fn.cluster.SCT(
  data = data,
  batch = "batch",
  vars.to.regress = c("G2M.Score", "S.Score", "percent.mt")
)

# -----------------------------
# 5) Remove ambient-like / low-quality clusters and re-run
# -----------------------------
data <- subset(
  data,
  subset = !(as.character(SCT_snn_res.0.2) %in% omit_clusters_res02)
)

data <- fn.cluster.SCT(
  data = data,
  batch = "batch",
  vars.to.regress = c("G2M.Score", "S.Score", "percent.mt")
)

# -----------------------------
# 6) Subclustering 1
# -----------------------------
Idents(data) <- data$SCT_snn_res.0.2
data.sub1 <- subset(data, idents = refine_clusters_1)

data.sub1 <- fn.cluster.SCT(
  data = data.sub1,
  batch = "batch",
  vars.to.regress = c("G2M.Score", "S.Score", "percent.mt")
)

map_sub1 <- c(
  # "0" = "Neutrophil",
  # "1" = "Omit"
)

if (length(map_sub1) == 0) {
  stop("map_sub1 is empty. Define subcluster->label mapping for subclustering 1 before running.")
}

data.sub1$cell_fine_sub1 <- unname(map_sub1[as.character(data.sub1$SCT_snn_res.0.2)])

data$cell_fine_sub1 <- NA_character_
data$cell_fine_sub1[Cells(data.sub1)] <- data.sub1$cell_fine_sub1

# -----------------------------
# 7) Subclustering 2
# -----------------------------
Idents(data) <- data$SCT_snn_res.0.2
data.sub2 <- subset(data, idents = refine_clusters_2)

data.sub2 <- fn.cluster.SCT(
  data = data.sub2,
  batch = "batch",
  vars.to.regress = c("G2M.Score", "S.Score", "percent.mt")
)

map_sub2 <- c(
  # "0" = "T.NKT"
)

if (length(map_sub2) == 0) {
  stop("map_sub2 is empty. Define subcluster->label mapping for subclustering 2 before running.")
}

data.sub2$cell_fine_sub2 <- unname(map_sub2[as.character(data.sub2$SCT_snn_res.0.2)])

data$cell_fine_sub2 <- NA_character_
data$cell_fine_sub2[Cells(data.sub2)] <- data.sub2$cell_fine_sub2

# -----------------------------
# 8) Final labels
# -----------------------------
Idents(data) <- data$SCT_snn_res.0.2
main_map_fine <- c(
  # "0" = "TAM.M2"
)

if (length(main_map_fine) == 0) {
  stop("main_map_fine is empty. Define main cluster->cell_fine mapping before running.")
}

data$cell_fine <- unname(main_map_fine[as.character(data$SCT_snn_res.0.2)])

ix1 <- which(!is.na(data$cell_fine_sub1))
data$cell_fine[ix1] <- data$cell_fine_sub1[ix1]
ix2 <- which(!is.na(data$cell_fine_sub2))
data$cell_fine[ix2] <- data$cell_fine_sub2[ix2]

data$cell_middle <- dplyr::recode(
  data$cell_fine,
  "Neutrophil"          = "Neutrophil",
  "Neutrophil.immature" = "Neutrophil",
  "Macrophage"          = "Macrophage",
  "TAM.M1"              = "Macrophage",
  "TAM.M2"              = "Macrophage",
  "Macrophage.resident" = "Macrophage",
  "Macrophage.ND"       = "Macrophage",
  "DC"                  = "DC",
  "Fibroblast"          = "Fibroblast",
  "Endothelial"         = "Endothelial",
  "SMC"                 = "SMC",
  "T.NKT"               = "T.NKT",
  "B"                   = "B",
  .default              = "Undetermined"
)

data$cell_coarse <- dplyr::recode(
  data$cell_middle,
  "Neutrophil"  = "Myeloid",
  "Macrophage"  = "Myeloid",
  "DC"          = "Myeloid",
  "T.NKT"       = "Lymphoid",
  "B"           = "Lymphoid",
  "Fibroblast"  = "Stromal",
  "Endothelial" = "Stromal",
  "SMC"         = "Stromal",
  .default      = "Undetermined"
)

# -----------------------------
# 9) Save
# -----------------------------
saveRDS(data, file.path(OUT_DIR, "PDX_mm_Integrated.rds"))
