#!/usr/bin/env Rscript
############################################################
# NicheNet Ligand Activity Analysis
#
# Input:
#   1) scRNAseq_PDX_NeuDep.rds  (integrated Seurat object)
#   2) DEGs_Meso_Isotype_vs_Ly6G.tsv
#
# Sender:
#   Neutrophils (Isotype condition)
#
# Receiver:
#   Tumor cells (Isotype vs Ly6G contrast)
#
# Output:
#   1) "NicheNetLigand_Priority_Specificity.tsv" # Table S7
#
# NicheNet version:
#   nichenetr 1.1.1
############################################################

suppressPackageStartupMessages({
  library(Seurat)        # v4
  library(nichenetr)     # 1.1.1
  library(dplyr)
  library(readr)
  library(tibble)
})

set.seed(123)

# ==========================================================
# 0) Paths (generalized)
# ==========================================================
BASE_DIR <- "/PATH/TO/PROJECT"

SEURAT_RDS <- file.path(BASE_DIR, "data", "scRNAseq_PDX_NeuDep.rds")
DEG_TSV    <- file.path(BASE_DIR, "data", "DEGs_Meso_Isotype_vs_Ly6G.tsv")

PRIOR_DIR  <- file.path(BASE_DIR, "data", "nichenet_priors")
LIGAND_TARGET_RDS     <- file.path(PRIOR_DIR, "ligand_target_matrix.rds")
WEIGHTED_NETWORKS_RDS <- file.path(PRIOR_DIR, "weighted_networks.rds")

OUT_DIR <- file.path(BASE_DIR, "output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
out_tsv = file.path(OUT_DIR, "NicheNetLigand_Priority_Specificity.tsv")

# ==========================================================
# 1) Load data and NicheNet priors
# ==========================================================
d <- readRDS(SEURAT_RDS)

ligand_target_matrix <- readRDS(LIGAND_TARGET_RDS)
weighted_networks    <- readRDS(WEIGHTED_NETWORKS_RDS)

lr_network <- readRDS(
  url("https://zenodo.org/record/3260758/files/lr_network.rds")
)

# ==========================================================
# 2) Define sender population (Neutrophils; Isotype)
# ==========================================================
d.neu <- subset(
  d,
  subset = condition == "Isotype" &
           cell_middle == "Neutrophil"
)

cell_rate <- 0.20

expr_sender <- GetAssayData(d.neu, assay = "RNA", slot = "data")
expr_prop_sender <- apply(expr_sender > 0, 1, mean)

expressed_genes_sender <- names(expr_prop_sender)[
  expr_prop_sender >= cell_rate
]

# ==========================================================
# 3) Define receiver population (Tumor cells)
# ==========================================================
d.meso <- subset(
  d,
  subset = cell_middle == "Tumor"
)

expr_receiver <- GetAssayData(d.meso, assay = "RNA", slot = "counts")
expr_prop_receiver <- apply(expr_receiver > 0, 1, mean)

expressed_genes_receiver <- names(expr_prop_receiver)[
  expr_prop_receiver >= cell_rate
]

background_expressed_genes <- intersect(
  expressed_genes_receiver,
  rownames(ligand_target_matrix)
)

# ==========================================================
# 4) Define gene set of interest (DE genes)
# ==========================================================
deg <- read_tsv(DEG_TSV, show_col_types = FALSE)

deg_filtered <- deg %>%
  filter(p_val_adj < 0.05,
         abs(avg_log2FC) > 0.5)

geneset_oi <- deg_filtered$Gene %>%
  intersect(rownames(ligand_target_matrix)) %>%
  intersect(background_expressed_genes)

# ==========================================================
# 5) Define potential ligands and receptors
# ==========================================================
ligands   <- unique(lr_network$from)
receptors <- unique(lr_network$to)

expressed_ligands  <- intersect(ligands, expressed_genes_sender)
expressed_receptors <- intersect(receptors, expressed_genes_receiver)

lr_network_expressed <- lr_network %>%
  filter(from %in% expressed_ligands,
         to   %in% expressed_receptors)

potential_ligands <- unique(lr_network_expressed$from)

# ==========================================================
# 6) Ligand activity analysis (NicheNet)
# ==========================================================
ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
)

ligand_activities <- ligand_activities %>%
  arrange(desc(pearson))

write_tsv(
  ligand_activities,
  file.path(OUT_DIR, "NicheNet_ligand_activity.tsv")
)

# ==========================================================
# 7) Ligand prioritization (cell-type specificity Z-score)
# ==========================================================

Idents(d) <- d$cell_middle

ligands_tested <- ligand_activities$test_ligand

# -----------------------------
# 7.1 Fraction of expressing cells
# -----------------------------
expr_all <- GetAssayData(d, assay = "RNA", slot = "data")
expr_all_sub <- expr_all[ligands_tested, , drop = FALSE]

cell_expr_fraction <- apply(expr_all_sub > 0, 1, function(x) {
  tapply(x, Idents(d), mean)
}) %>% t()

# -----------------------------
# 7.2 Non-negative average expression
#     (library-size normalized counts)
# -----------------------------
NonNegativeAverageExpression <- function(seurat_obj,
                                         group_by = "cell_middle",
                                         scale.factor = 10000) {

  counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
  groups <- seurat_obj[[group_by]][,1]

  count_aggr <- apply(counts, 1, function(gene_counts) {
    tapply(gene_counts, groups, sum)
  }) %>% t()

  libsize_group <- colSums(count_aggr)

  scaled <- sweep(count_aggr, 2, libsize_group, "/") * scale.factor
  log1p(scaled)
}

log_avg <- NonNegativeAverageExpression(d)
log_avg_sub <- log_avg[ligands_tested, , drop = FALSE]

# -----------------------------
# 7.3 Ligand priority score
# -----------------------------
ligand_priority <- log_avg_sub * cell_expr_fraction

# -----------------------------
# 7.4 Z-score across cell types
# -----------------------------
ligand_priority_z <- t(apply(ligand_priority, 1, function(x) {
  if (sd(x) == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}))

colnames(ligand_priority_z) <- paste0(colnames(ligand_priority_z), ".z.score")

# -----------------------------
# 7.5 Export
# -----------------------------
write_tsv(
  as.data.frame(ligand_priority_z) %>%
    rownames_to_column("Ligand"),
  file.path(OUT_DIR, out_tsv)
)
