#!/usr/bin/env Rscript

############################################################
# 03_prepare_stage_histology_collapse.R
#
# Purpose:
#   Prepare PanCan survival table with:
#     - HA.Strata defined within each acronym (Q4 vs Q1-3)
#     - pathologic_stage collapsed using a lookup TSV
#     - histological_type collapsed using a lookup TSV
#     - remove sparse categories (n < 5) within each acronym
#
# Inputs (EDIT BEFORE RUN):
#   - TSV created in earlier step 
#     Must include: acronym, days_to_censor, status, HA
#     Optionally: age_at_initial_pathologic_diagnosis, gender,
#                pathologic_stage, histological_type
#   - data/lookup_stage_collapse.tsv
#   - data/lookup_hist_collapse.tsv
#
# Outputs:
#   - output/pancan_surv_table_collapsed.tsv
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
})

# -----------------------------
# Paths (EDIT BEFORE RUN)
# -----------------------------
BASE_DIR <- "/PATH/TO/TCGA"
DATA_DIR <- file.path(BASE_DIR, "data")

# Input RDS (created previously)
IN_TSV <- file.path(OUT_DIR, "Immune10_deconvolv_TCGA.PanCan.annotated.tsv")

# Lookup TSVs
IN_STAGE_LOOKUP <- file.path(DATA_DIR, "lookup_stage_collapse.tsv")
IN_HIST_LOOKUP  <- file.path(DATA_DIR, "lookup_hist_collapse.tsv")

# Outputs
OUT_RDS <- file.path(OUT_DIR, "pancan_surv_table_collapsed.rds")
OUT_TSV <- file.path(OUT_DIR, "pancan_surv_table_collapsed.tsv")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1) Load
# -----------------------------
df <- read_tsv(IN_TSV)
stopifnot(all(c("acronym", "days_to_censor", "status", "HA") %in% colnames(df)))

# -----------------------------
# 2) HA Strata (Q4 vs Q1-3) within each acronym
# -----------------------------
df <- df %>%
  group_by(acronym) %>%
  mutate(
    HA.Strata = ifelse(HA > quantile(HA, probs = 0.75, na.rm = TRUE), "High", "Low"),
    HA.Strata = factor(HA.Strata, levels = c("Low", "High"))
  ) %>%
  ungroup()

# -----------------------------
# 3) Basic cleanup + rename
# -----------------------------
if ("age_at_initial_pathologic_diagnosis" %in% colnames(df)) {
  df <- df %>% rename(Age = age_at_initial_pathologic_diagnosis)
}

df <- df %>%
  mutate(
    days_to_censor = as.numeric(days_to_censor),
    status = as.numeric(status)
  )

if ("Age" %in% colnames(df)) df <- df %>% mutate(Age = as.numeric(Age))
if ("gender" %in% colnames(df)) df <- df %>% mutate(gender = as.character(gender))
if ("pathologic_stage" %in% colnames(df)) df <- df %>% mutate(pathologic_stage = as.character(pathologic_stage))
if ("histological_type" %in% colnames(df)) df <- df %>% mutate(histological_type = as.character(histological_type))

# -----------------------------
# 4) Collapse pathologic_stage using lookup TSV
# -----------------------------
if ("pathologic_stage" %in% colnames(df)) {
  lookup_stage <- read_tsv(IN_STAGE_LOOKUP, show_col_types = FALSE) %>%
    select(stage_original, stage_collapsed)

  df <- df %>%
    left_join(lookup_stage, by = c("pathologic_stage" = "stage_original")) %>%
    mutate(
      pathologic_stage = ifelse(!is.na(stage_collapsed), stage_collapsed, pathologic_stage)
    ) %>%
    select(-stage_collapsed)

  # cancer-specific overrides (kept minimal; edit if needed)
  df$pathologic_stage[df$acronym %in% c("GBM", "LGG")] <- "[Not Applicable]"
  df$pathologic_stage[df$acronym == "CESC"] <- "NA"
}

# -----------------------------
# 5) Collapse histological_type using lookup TSV
# -----------------------------
if ("histological_type" %in% colnames(df)) {
  lookup_hist <- read_tsv(IN_HIST_LOOKUP, show_col_types = FALSE) %>%
    select(acronym, histological_type, new_hist)

  df <- df %>%
    left_join(lookup_hist, by = c("acronym", "histological_type")) %>%
    mutate(
      histological_type = ifelse(!is.na(new_hist), new_hist, histological_type)
    ) %>%
    select(-new_hist)
}

# -----------------------------
# 6) Filters for downstream multiv analyses
# -----------------------------
df <- df %>%
  filter(!is.na(days_to_censor), !is.na(status), !is.na(HA.Strata))

# Optional: drop not-available stage
if ("pathologic_stage" %in% colnames(df)) {
  df <- df %>% filter(pathologic_stage != "[Not Available]")
}

# Remove sparse histology within each acronym (n < 5)
if ("histological_type" %in% colnames(df)) {
  df <- df %>%
    group_by(acronym, histological_type) %>%
    filter(n() >= 5) %>%
    ungroup()
}

# -----------------------------
# 7) Save
# -----------------------------
write_tsv(df, OUT_TSV)

message("Saved:")
message("  ", OUT_TSV)

