#!/usr/bin/env Rscript

# ============================================================
# PDX human tumor (hg) integration: Isotype vs Ly6G
# Seurat v4.3.0.1
# ============================================================

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

# -----------------------------
# 0) Resolve repo paths (stable)
# -----------------------------
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
SCRIPT_DIR <- normalizePath(dirname(script_path))
REPO_ROOT  <- normalizePath(file.path(SCRIPT_DIR, ".."))  # scRNA_processing/

source(file.path(REPO_ROOT, "R", "fn.cluster.SCT.R"))

data("cc.genes", package = "Seurat")

# -----------------------------
# 1) Inputs (generalized)
# -----------------------------
BASE_DIR   <- "/PATH/TO/SEURAT_OBJECTS/"
OUT_DIR    <- file.path(REPO_ROOT, "output", "PDX_integration_hg")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

HA_DEG_TSV <- "/PATH/TO/DEG_HA_vs_Control.tsv" # Table S8

sampleIDs <- c("PDX.M_Isotype", "PDX.M_Ly6G")

# -----------------------------
# 2) Load and merge
# -----------------------------
obj_list <- lapply(sampleIDs, function(sid) {
  readRDS(file.path(BASE_DIR, paste0(sid, "_hg.rds")))
})

data <- merge(obj_list[[1]], y = obj_list[[2]])

data$condition <- ifelse(grepl("Isotype", data$batch), "Isotype", "Ly6G")
Idents(data) <- data$condition

# (以下、解析ロジックはそのまま)
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

data <- CellCycleScoring(
  data,
  s.features   = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes
)

data <- fn.cluster.SCT(
  data = data,
  batch = "condition",
  vars.to.regress = c("G2M.Score", "S.Score", "percent.mt")
)

data <- subset(data, subset = SCT_snn_res.1.2 != 7)

data <- fn.cluster.SCT(
  data = data,
  batch = "condition",
  vars.to.regress = c("G2M.Score", "S.Score", "percent.mt")
)

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

deg_HA <- read_tsv(HA_DEG_TSV)

genes_HA_UP <- deg_HA %>%
  filter(logFC > 0 & adj.P.Val < 0.05) %>%
  pull(Gene)

data <- AddModuleScore(
  object   = data,
  features = list(genes_HA_UP),
  name     = "HA.UP.Score"
)

saveRDS(data, file.path(OUT_DIR, "PDX_hg_Integrated.rds"))
