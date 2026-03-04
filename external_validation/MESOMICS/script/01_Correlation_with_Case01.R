#!/usr/bin/env Rscript

############################################################
# 01_Correlation_with_Case01.R
#
# Purpose:
#   Correlate Case01 tumor pseudo-bulk (scRNA-seq) with MESOMICS bulk RNA-seq.
#
# Inputs (generalized):
#   - Case01 Seurat object (integrated / annotated): scRNAseq_Case01.rds
#   - MESOMICS gene-level counts table (GeneSymbol): MESOMICS_Count_GeneSymbol.tsv
#   - MESOMICS archetypes table: TableS28_Archetypes.tsv
#   - MESOMICS sample overview: TableS2-3_SamplesOverview.tsv
#   - BAP1 mutation/CNV table:
#   - MESOMICS CNV summary (TERT): TableS36_AMP.DEL.genes.tsv
#
# Outputs:
#   output/Corr_with_Case01.psbulk.tsv
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(quantiseqr)
})

# -----------------------------
# Paths (EDIT BEFORE RUN)
# -----------------------------
BASE_DIR <- "/PATH/TO/external_validation/MESOMICS"

IN_SEURAT_RDS       <- "/PATH/TO/scRNAseq_Case01.rds"
IN_MESOMICS_COUNTS  <- file.path(BASE_DIR, "data", "MESOMICS_Count_GeneSymbol.tsv")
IN_ARCHETYPES       <- file.path(BASE_DIR, "data", "TableS28_Archetypes.tsv")
IN_SAMPLES_OVERVIEW <- file.path(BASE_DIR, "data", "TableS2-3_SamplesOverview.tsv")
IN_BAP1_MUT         <- file.path(BASE_DIR, "data", "BAP1.SNV.CNV_c.tsv")
IN_TERT_CNV         <- file.path(BASE_DIR, "data", "TableS36_AMP.DEL.genes.tsv")

OUT_DIR <- file.path(BASE_DIR, "output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_TSV <- file.path(OUT_DIR, "Corr_with_Case01.psbulk.tsv")

# -----------------------------
# Helpers
# -----------------------------
fn.perMillion <- function(df) {
  as.data.frame(apply(df, 2, function(x) x / sum(as.numeric(x)) * 1e6))
}

# -----------------------------
# Load data
# -----------------------------
d <- readRDS(IN_SEURAT_RDS)

d.mesomics.tpm <- read_tsv(IN_MESOMICS_COUNTS, show_col_types = FALSE) %>%
  column_to_rownames("Symbol") %>%
  fn.perMillion()

d.arc  <- read_tsv(IN_ARCHETYPES, show_col_types = FALSE)
d.anno <- read_tsv(IN_SAMPLES_OVERVIEW, show_col_types = FALSE)

d.mut <- read_tsv(IN_BAP1_MUT, show_col_types = FALSE) %>%
  replace_na(replace = list(Variant_Classification = "No Mutation")) %>%
  mutate(
    BAP1_status = ifelse(
      (CNV_BAP1 %in% c("Homozygous deletion", "Heterozygous deletion") |
         Variant_Classification != "No Mutation"),
      "altered", "unaltered"
    )
  )

# C/V score (CLDN15 / VIM)
cv <- {
  C <- d.mesomics.tpm["CLDN15", ]
  V <- d.mesomics.tpm["VIM", ]
  (log2(C + 1) / log2(V + 1)) %>%
    t() %>% as.data.frame() %>%
    rownames_to_column("ID") %>%
    rename(cv = CLDN15)
}

# TERT status
d.tert <- read_tsv(IN_TERT_CNV, show_col_types = FALSE) %>%
  select(Cohort, TERT) %>%
  mutate(
    TERT = TERT %>%
      gsub("Amplification", "Amp", .) %>%
      gsub("Heterozygous deletion", "Het Del", .) %>%
      gsub("No CN change", "Neutral", .),
    TERT = factor(TERT, levels = c("Amp", "Neutral", "nLOH", "Het Del"))
  )

# -----------------------------
# Annotation merge (inner joins; may reduce sample count)
# -----------------------------
ann <- merge(d.mut, d.anno, by.x = "Tumor_Sample_Barcode", by.y = "Sample")
ann <- merge(d.arc, ann, by.x = "ID", by.y = "Tumor_Sample_Barcode")
ann <- merge(ann, d.tert, by.x = "ID", by.y = "Cohort")
ann <- merge(cv, ann, by = "ID")
ann <- ann[!duplicated(ann$ID), ]

ann <- ann %>%
  mutate(
    Type = Type %>% gsub("MME", "PM.E", .) %>% gsub("MMB", "PM.B", .) %>% gsub("MMS", "PM.S", .),
    Type = factor(Type, levels = c("PM.E", "PM.B", "PM.S")),
    IHC.BAP1 = IHC.BAP1 %>% gsub("YES in MMS/NO in MME", "Loss", .),
    IHC.BAP1 = IHC.BAP1 %>% gsub("YES", "Retain", .) %>% gsub("NO", "Loss", .),
    IHC.BAP1 = factor(IHC.BAP1, levels = c("Retain", "Loss"))
  )

# keep common samples
comm <- intersect(ann$ID, colnames(d.mesomics.tpm))
ann <- ann[match(comm, ann$ID), , drop = FALSE]
d.mesomics.tpm <- d.mesomics.tpm[, match(comm, colnames(d.mesomics.tpm)), drop = FALSE]
stopifnot(all(ann$ID == colnames(d.mesomics.tpm)))

# -----------------------------
# quanTIseq (optional metadata)
# -----------------------------
ti_mesomics <- quantiseqr::run_quantiseq(
  expression_data = d.mesomics.tpm,
  signature_matrix = "TIL10",
  is_arraydata = FALSE,
  is_tumordata = TRUE,
  scale_mRNA = TRUE
)

# -----------------------------
# Case01 tumor pseudo-bulk
# -----------------------------
d.tum <- subset(d, subset = cell_coarse == "Mesothelial")
d.tum$cell_coarse <- factor(d.tum$cell_coarse, levels = "Mesothelial")

d.bulk <- PseudobulkExpression(d.tum, return.seurat = FALSE, group.by = "cell_coarse")

d.case01.tpm <- fn.perMillion(d.bulk$RNA) %>%
  rownames_to_column("Symbol") %>%
  rename(Case01_Tumor = all)

# -----------------------------
# Correlation across genes (log2(per-million + 1))
# -----------------------------
dm <- merge(
  d.case01.tpm,
  d.mesomics.tpm %>% rownames_to_column("Symbol"),
  by = "Symbol"
) %>% column_to_rownames("Symbol")

dm <- dm[dm$Case01_Tumor != 0, , drop = FALSE]  # dropout filter
dm <- fn.perMillion(dm)
dml <- log2(dm + 1)

cor_vec <- apply(dml, 2, function(x) round(cor(dml$Case01_Tumor, x), 3))
cors <- tibble(sample = names(cor_vec), Corr_with_Case01.psbulk = as.numeric(cor_vec)) %>%
  filter(sample != "Case01_Tumor") %>%
  left_join(ann, by = c("sample" = "ID")) %>%
  mutate(Tertile = ntile(Corr_with_Case01.psbulk, 3),
         Tertile = factor(ifelse(Tertile == 1, "Low", ifelse(Tertile == 2, "Moderate", "High")),
                          levels = c("High", "Moderate", "Low")))

# add quanTIseq (matched by sample colnames)
ti_df <- as.data.frame(ti_mesomics) %>% rownames_to_column("sample")
cors <- cors %>% left_join(ti_df, by = "sample")

# export
cors %>%
  rename(`MOFA.MESOMICS.Cell.division` = `MOFA.MESOMICS.Cell division`) %>%
  write_tsv(OUT_TSV)

message("Wrote: ", OUT_TSV)

