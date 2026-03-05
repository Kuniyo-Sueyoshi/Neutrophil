#!/usr/bin/env Rscript

############################################################
# 03_PDX_sample_integration_hg.R
#
# Purpose:
#   PDX human tumor (hg) integration: Isotype vs Ly6G.
#   - Merge hg compartments
#   - QC filtering, cell cycle scoring
#   - SCT-based integration/clustering
#   - DEG (Isotype vs Ly6G)
#   - GO enrichment (Isotype-up genes)
#   - Module scores (migration/chemotaxis and HA-response genes)
#
# Inputs (repo-local):
#   output/PDX_split/PDX.M_Isotype/PDX.M_Isotype_hg.rds
#   output/PDX_split/PDX.M_Ly6G/PDX.M_Ly6G_hg.rds
#
# External input (EDIT):
#   HA_DEG_TSV : DEG table for histamine-response genes (e.g., Table S8)
#
# Outputs (repo-local):
#   output/PDX_integration_hg/DEGs_Meso_Isotype_vs_Ly6G.tsv
#   output/PDX_integration_hg/GO_BP_enrichment_Isotype_UP.tsv
#   output/PDX_integration_hg/PDX_hg_Integrated.rds
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(tibble)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

set.seed(123)

# ==========================
# Paths (EDIT BEFORE RUN)
# ==========================
BASE_DIR <- "/PATH/TO/scRNA_processing"

OUT_BASE <- file.path(BASE_DIR, "output")
R_DIR    <- file.path(BASE_DIR, "R")

# External input (EDIT)
HA_DEG_TSV <- "/PATH/TO/DEG_HA_vs_Control.tsv"  # e.g., Table S8

# Output
OUT_DIR <- file.path(OUT_BASE, "PDX_integration_hg")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ==========================
# Inputs (repo-local)
# ==========================
sampleIDs <- c("PDX.M_Isotype", "PDX.M_Ly6G")

IN_HG_RDS <- setNames(
  file.path(OUT_BASE, "PDX_split", sampleIDs, paste0(sampleIDs, "_hg.rds")),
  sampleIDs
)

# ==========================
# Source helpers
# ==========================
source(file.path(R_DIR, "fn.cluster.SCT.R"))
data("cc.genes", package = "Seurat")

# ============================================================
# 1) Load and merge
# ============================================================
obj_list <- lapply(sampleIDs, function(sid) readRDS(IN_HG_RDS[[sid]]))
data <- merge(obj_list[[1]], y = obj_list[[2]])

# condition label from batch (kept as original logic)
data$condition <- ifelse(grepl("Isotype", data$batch), "Isotype", "Ly6G")
Idents(data) <- data$condition

# ============================================================
# 2) QC filtering (same thresholds as original)
# ============================================================
UMI.cutoff     <- 500
Feature.cutoff <- 200
Mt.cutoff      <- 15
Hb.cutoff      <- 5
ribo.cutoff    <- 5

data <- subset(
  data,
  subset =
    nCount_RNA    > UMI.cutoff &
    nFeature_RNA  > Feature.cutoff &
    percent.mt    < Mt.cutoff &
    percent.hb    < Hb.cutoff &
    percent.ribo  > ribo.cutoff
)

# ============================================================
# 3) Cell cycle scoring
# ============================================================
data <- CellCycleScoring(
  data,
  s.features   = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes
)

# ============================================================
# 4) SCT clustering (Isotype vs Ly6G integration)
# ============================================================
data <- fn.cluster.SCT(
  data = data,
  batch = "condition",
  vars.to.regress = c("G2M.Score", "S.Score", "percent.mt")
)

# remove cluster 7 (as in original)
data <- subset(data, subset = SCT_snn_res.1.2 != 7)

data <- fn.cluster.SCT(
  data = data,
  batch = "condition",
  vars.to.regress = c("G2M.Score", "S.Score", "percent.mt")
)

# ============================================================
# 5) DEG: Isotype vs Ly6G
# ============================================================
Idents(data) <- data$condition

deg <- FindMarkers(
  object = data,
  ident.1 = "Isotype",
  ident.2 = "Ly6G",
  min.pct = 0.1,
  logfc.threshold = 0,
  pseudocount.use = 0.1
) %>%
  rownames_to_column("Gene") %>%
  filter(!grepl("^(MT-|mt-)|^Gm[0-9]+", Gene))

write_tsv(deg, file.path(OUT_DIR, "DEGs_Meso_Isotype_vs_Ly6G.tsv"))

# ============================================================
# 6) GO enrichment: Isotype-up genes
# ============================================================
gene_up <- deg %>%
  filter(avg_log2FC > 0.5 & p_val_adj < 0.05) %>%
  pull(Gene)

ego <- enrichGO(
  gene      = gene_up,
  OrgDb     = org.Hs.eg.db,
  keyType   = "SYMBOL",
  ont       = "BP",
  minGSSize = 30,
  maxGSSize = 300
)

ego <- simplify(ego)
ego_res <- ego@result %>% arrange(p.adjust)

write_tsv(ego_res, file.path(OUT_DIR, "GO_BP_enrichment_Isotype_UP.tsv"))

# ============================================================
# 7) Module scores
# ============================================================
mig_genes <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = "GO:0010634",
  keytype = "GOALL",
  columns = "SYMBOL"
) %>% distinct(SYMBOL)

chemo_genes <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = "GO:1990266",
  keytype = "GOALL",
  columns = "SYMBOL"
) %>% distinct(SYMBOL)

data <- AddModuleScore(
  object   = data,
  features = list(mig_genes$SYMBOL, chemo_genes$SYMBOL),
  name     = c("Migration.Score", "Chemotaxis.Score")
)

deg_HA <- read_tsv(HA_DEG_TSV, show_col_types = FALSE)
genes_HA_UP <- deg_HA %>%
  filter(logFC > 0 & adj.P.Val < 0.05) %>%
  pull(Gene)

data <- AddModuleScore(
  object   = data,
  features = list(genes_HA_UP),
  name     = "HA.UP.Score"
)

# ============================================================
# 8) Save integrated object
# ============================================================
saveRDS(data, file.path(OUT_DIR, "PDX_hg_Integrated.rds"))

