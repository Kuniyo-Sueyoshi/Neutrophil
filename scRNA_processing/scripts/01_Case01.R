#####################################################
# 01_Case01_preprocess.R
# scRNA-seq preprocessing and analysis pipeline (Case01)
#
# Purpose:
# Preprocessing + analysis script for Case01:
# QC metrics, SoupX ambient RNA correction, QC filtering,
# doublet removal, cell cycle scoring (metadata only),
# clustering (LogNormalize), tumor ribosomal filtering,
# DEG, and pathway enrichment (fgsea).
#
# Inputs (generalized):
# 10x Genomics output directory containing:
# /PATH/TO/10X_OUTPUTS/Case01_T/outs/
# - filtered_feature_bc_matrix.h5
# - raw_feature_bc_matrix.h5 (required for SoupX)
#
# Outputs:
# - Processed Seurat object (RDS, scRNAseq_Case01.rds)
# - DEG table (TSV, Table S1) # - GSEA/fgsea results (TSV, Table S2)
##########################################################

############################################################
# 01_Case01_preprocess.R
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(SoupX)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(fgsea)
})

set.seed(123)

# -----------------------------
# User parameters (generalized)
# -----------------------------
HOME_DIR  <- "/PATH/TO/10X_OUTPUTS"
SAMPLE_ID <- "Case01_T"

# output: repo-local (scRNA_processing/output/<SAMPLE_ID>)
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
SCRIPT_DIR <- normalizePath(dirname(script_path))
REPO_ROOT  <- normalizePath(file.path(SCRIPT_DIR, ".."))  # scRNA_processing/
OUT_DIR    <- file.path(REPO_ROOT, "output", SAMPLE_ID)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

QC <- list(
  umi_min        = 500,
  feat_min       = 200,
  mt_max         = 15,
  hb_max         = 5,
  ribo_min_tumor = 5
)

MSIGDB_HALLMARK_GMT <- "/PATH/TO/MSIGDB/h.all.v2024.1.Hs.symbols.gmt"

# -----------------------------
# Helper functions
# -----------------------------
# source paths: stable regardless of working directory
source(file.path(REPO_ROOT, "R", "fn.cluster.log.R"))
source(file.path(REPO_ROOT, "R", "fn.QCmetrics.R"))

load_filtered_h5 <- function(home_dir, sample_id) {
  f <- file.path(home_dir, sample_id, "outs", "filtered_feature_bc_matrix.h5")
  mat <- Read10X_h5(f)
  CreateSeuratObject(counts = mat, project = sample_id)
}

run_soupx <- function(home_dir, sample_id) {
  dir_10x <- file.path(home_dir, sample_id, "outs")
  sc <- load10X(dir_10x)
  sc <- autoEstCont(sc, forceAccept = TRUE)
  adj <- adjustCounts(sc, roundToInt = TRUE)
  list(sc = sc, counts = adj)
}

run_doublet_filter <- function(seurat_obj) {
  set.seed(123)
  sce <- scDblFinder(as.SingleCellExperiment(seurat_obj))
  obj <- as.Seurat(sce, data = NULL)
  subset(obj, subset = scDblFinder.class == "singlet")
}

read_gmt <- function(gmt_file) {
  raw <- scan(gmt_file, what = "", sep = "\n", quiet = TRUE)
  spl <- strsplit(raw, "\t")
  pathways <- lapply(spl, function(x) x[-c(1, 2)])
  names(pathways) <- vapply(spl, `[`, character(1), 1)
  pathways
}

# ============================================================
# 1) Load filtered counts and compute QC metrics (pre-SoupX)
# ============================================================
obj0 <- load_filtered_h5(HOME_DIR, SAMPLE_ID)

# NOTE: keep the function call style consistent.
# If add_qc_metrics() expects (obj, sample_id, cap=...), use the same signature everywhere.
obj0 <- add_qc_metrics(obj0, cap = FALSE)

# BUGFIX: sampleID -> SAMPLE_ID
obj0$batch <- SAMPLE_ID

# ============================================================
# 2) SoupX ambient RNA correction (evaluated)
# ============================================================
soupx <- run_soupx(HOME_DIR, SAMPLE_ID)

obj <- CreateSeuratObject(
  counts = soupx$counts,
  project = SAMPLE_ID,
  min.cells = 1,
  min.features = 1
)

# keep consistent signature
obj <- add_qc_metrics(obj, cap = FALSE)

# ============================================================
# 3) Global QC filtering
# ============================================================
obj <- subset(
  obj,
  subset =
    nCount_RNA   > QC$umi_min &
    nFeature_RNA > QC$feat_min &
    percent.mt   < QC$mt_max &
    percent.hb   < QC$hb_max
)

# ============================================================
# 4) Doublet detection/removal
# ============================================================
obj <- run_doublet_filter(obj)
obj$batch <- SAMPLE_ID

# ============================================================
# 5) Cell cycle scoring
# ============================================================
s.genes   <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes

obj <- Seurat::CellCycleScoring(
  object = obj,
  s.features = s.genes,
  g2m.features = g2m.genes,
  verbose = FALSE
)

# ============================================================
# 6) Initial clustering
# ============================================================
obj <- fn.cluster.log(
  data = obj,
  vars.to.regress = c("percent.mt", "nCount_RNA"),
  resolutions = c(0.2, 0.4)
)

# ============================================================
# 8) Tumor-specific ribosomal RNA filtering
# ============================================================
obj <- subset(
  obj,
  subset =
    (cell_coarse != "Mesothelial") |
    (cell_coarse == "Mesothelial" & percent.ribo > QC$ribo_min_tumor)
)

# ============================================================
# 9) Re-cluster after tumor ribosomal filter
# ============================================================
obj <- fn.cluster.log(
  data = obj,
  vars.to.regress = c("percent.mt", "nCount_RNA"),
  resolutions = c(0.2, 0.4, 0.8)
)

# ============================================================
# 11) Differential expression (Tumor.M vs Tumor.K)
# ============================================================
obj_tumor <- subset(obj, subset = !is.na(TumorSubtype))
Idents(obj_tumor) <- obj_tumor$TumorSubtype

deg <- FindMarkers(
  object = obj_tumor,
  ident.1 = "Tumor.M",
  ident.2 = "Tumor.K",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0.01
) |>
  rownames_to_column("Gene")

write.table(
  deg,
  file = file.path(OUT_DIR, "DEG_TumorM_vs_TumorK.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ============================================================
# 12) Pathway enrichment (fgsea; Hallmark)
# ============================================================
stats <- deg$avg_log2FC
names(stats) <- deg$Gene
stats <- sort(stats, decreasing = FALSE)

pathways <- read_gmt(MSIGDB_HALLMARK_GMT)

gsea_res <- fgsea(
  pathways = pathways,
  stats    = stats,
  minSize  = 10,
  maxSize  = Inf
)

write.table(
  gsea_res,
  file = file.path(OUT_DIR, "GSEA_Hallmark.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ============================================================
# 13) Save processed object
# ============================================================
saveRDS(obj, file = file.path(OUT_DIR, paste0(SAMPLE_ID, "_processed.rds")))
