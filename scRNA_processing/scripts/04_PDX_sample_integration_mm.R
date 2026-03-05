#!/usr/bin/env Rscript

############################################################
# 04_PDX_sample_integration_mm.R
#
# Purpose:
#   PDX mouse compartment integration (PDX.M.Isotype vs PDX.M.Ly6G)
#
# Steps:
#   1) Load + merge mm Seurat objects
#   2) QC metrics annotation
#   3) Cell-cycle scoring (human cc.genes → mouse orthologs via gorth)
#   4) SCTransform + Harmony integration + clustering
#   5) Remove ambient-like clusters and re-run
#   6) Subclustering for ambiguous populations
#   7) Final cell-type labels
#   8) Save integrated object
#
# Inputs (repo-local):
#   output/PDX_split/PDX.M_Isotype/PDX.M_Isotype_mm.rds
#   output/PDX_split/PDX.M_Ly6G/PDX.M_Ly6G_mm.rds
#
# Output:
#   output/PDX_integration_mm/PDX_mm_Integrated.rds
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(gprofiler2)
})

set.seed(123)

# ==========================
# Paths (EDIT BEFORE RUN)
# ==========================
BASE_DIR <- "/PATH/TO/scRNA_processing"

OUT_BASE <- file.path(BASE_DIR, "output")
R_DIR    <- file.path(BASE_DIR, "R")

OUT_DIR <- file.path(OUT_BASE, "PDX_integration_mm")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ==========================
# Inputs
# ==========================
sampleIDs <- c("PDX.M_Isotype", "PDX.M_Ly6G")

in_files <- file.path(
  OUT_BASE,
  "PDX_split",
  sampleIDs,
  paste0(sampleIDs, "_mm.rds")
)

# ==========================
# Parameters
# ==========================
omit_clusters_res02 <- c("5", "13")
refine_clusters_1   <- c("1", "4")
refine_clusters_2   <- c("3", "12")

# ==========================
# Helper functions
# ==========================
source(file.path(R_DIR, "fn.cluster.SCT.R"))
source(file.path(R_DIR, "fn.QCmetrics.R"))

# -----------------------------
# Sanity checks
# -----------------------------
for (f in in_files) stopifnot(file.exists(f))

# ============================================================
# 1) Load and merge mm objects
# ============================================================
data <- merge(readRDS(in_files[1]), readRDS(in_files[2]))

if (!"batch" %in% colnames(data@meta.data)) {
  stop("meta.data$batch is missing in the input objects.")
}

data$condition <- ifelse(
  grepl("Isotype", data$batch, ignore.case = TRUE),
  "Isotype",
  "Ly6G"
)

# ============================================================
# 2) QC metrics annotation
# ============================================================
data <- add_qc_metrics(data, cap = TRUE)

# ============================================================
# 3) Cell-cycle scoring (human → mouse orthologs)
# ============================================================
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

# ============================================================
# 4) Integration (SCT + Harmony)
# ============================================================
data <- fn.cluster.SCT(
  data = data,
  batch = "batch",
  vars.to.regress = c("G2M.Score", "S.Score", "percent.mt")
)

# ============================================================
# 5) Remove ambient-like clusters
# ============================================================
data <- subset(
  data,
  subset = !(as.character(SCT_snn_res.0.2) %in% omit_clusters_res02)
)

data <- fn.cluster.SCT(
  data = data,
  batch = "batch",
  vars.to.regress = c("G2M.Score", "S.Score", "percent.mt")
)

# ============================================================
# 6) Subclustering 1
# ============================================================
Idents(data) <- data$SCT_snn_res.0.2
data.sub1 <- subset(data, idents = refine_clusters_1)

data.sub1 <- fn.cluster.SCT(
  data = data.sub1,
  batch = "batch",
  vars.to.regress = c("G2M.Score", "S.Score", "percent.mt")
)

map_sub1 <- c(
  # "0" = "Neutrophil"
)

if (length(map_sub1) == 0) {
  stop("Define subcluster→label mapping for subclustering 1.")
}

data.sub1$cell_fine_sub1 <- unname(map_sub1[as.character(data.sub1$SCT_snn_res.0.2)])

data$cell_fine_sub1 <- NA_character_
data$cell_fine_sub1[Cells(data.sub1)] <- data.sub1$cell_fine_sub1

# ============================================================
# 7) Subclustering 2
# ============================================================
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
  stop("Define subcluster→label mapping for subclustering 2.")
}

data.sub2$cell_fine_sub2 <- unname(map_sub2[as.character(data.sub2$SCT_snn_res.0.2)])

data$cell_fine_sub2 <- NA_character_
data$cell_fine_sub2[Cells(data.sub2)] <- data.sub2$cell_fine_sub2

# ============================================================
# 8) Final labels
# ============================================================
Idents(data) <- data$SCT_snn_res.0.2

main_map_fine <- c(
  # "0" = "TAM.M2"
)

if (length(main_map_fine) == 0) {
  stop("Define main cluster→cell_fine mapping.")
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

# ============================================================
# 9) Save
# ============================================================
saveRDS(data, file.path(OUT_DIR, "PDX_mm_Integrated.rds"))

