#!/usr/bin/env Rscript

############################################################
# 04_pancan_multivar.R
#
# Purpose:
#   Pan-cancer multivariable Cox analysis using the collapsed table:
#     - HA.Strata (Q4 vs Q1-3; within acronym)
#     - Age
#     - gender (if variable within acronym)
#     - pathologic_stage (collapsed; if available/variable within acronym)
#     - histological_type (collapsed; if available/variable within acronym)
#
# Input:
#   - output/pancan_surv_table_collapsed.tsv
#
# Outputs:
#   - output/pancan_multivar_cox_results.tsv  (all terms per acronym)
#   - output/pancan_multivar_cox_HA_only.tsv  (HA.Strata term per acronym)
#
# Notes:
#   - For each acronym, categorical covariates with only 1 level are excluded.
#   - Reference levels are set using data/lookup_hist_reference.tsv when possible.
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(broom)
})

# -----------------------------
# Paths (EDIT BEFORE RUN)
# -----------------------------
BASE_DIR <- "/PATH/TO/TCGA"
DATA_DIR <- file.path(BASE_DIR, "data")
OUT_DIR  <- file.path(BASE_DIR, "output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

IN_TSV <- file.path(OUT_DIR, "pancan_surv_table_collapsed.tsv")

IN_HIST_REF <- file.path(DATA_DIR, "lookup_hist_reference.tsv")  # acronym, hist_ref (optional but recommended)

OUT_ALL_TSV <- file.path(OUT_DIR, "pancan_multivar_cox_results.tsv")
OUT_HA_TSV  <- file.path(OUT_DIR, "pancan_multivar_cox_HA_only.tsv")

# -----------------------------
# 1) Load
# -----------------------------
df <- read_tsv(IN_TSV, show_col_types = FALSE)

stopifnot(all(c("acronym", "days_to_censor", "status", "HA.Strata") %in% colnames(df)))

# Standardize types
df <- df %>%
  mutate(
    acronym = as.character(acronym),
    days_to_censor = as.numeric(days_to_censor),
    status = as.numeric(status),
    HA.Strata = factor(HA.Strata, levels = c("Low", "High"))
  )

if ("Age" %in% colnames(df)) df <- df %>% mutate(Age = as.numeric(Age))
if ("gender" %in% colnames(df)) df <- df %>% mutate(gender = as.factor(gender))
if ("pathologic_stage" %in% colnames(df)) df <- df %>% mutate(pathologic_stage = as.factor(pathologic_stage))
if ("histological_type" %in% colnames(df)) df <- df %>% mutate(histological_type = as.factor(histological_type))

# Optional: reference table for histology
hist_ref <- NULL
if (file.exists(IN_HIST_REF)) {
  hist_ref <- read_tsv(IN_HIST_REF, show_col_types = FALSE) %>%
    select(acronym, hist_ref) %>%
    mutate(acronym = as.character(acronym), hist_ref = as.character(hist_ref))
}

# -----------------------------
# 2) Helper: run per-acronym Cox
# -----------------------------
run_cox_one_acronym <- function(d, hist_ref_tbl = NULL) {
  cancer <- unique(d$acronym)
  stopifnot(length(cancer) == 1)

  # candidate covariates
  base_vars <- c("HA.Strata")
  if ("Age" %in% colnames(d)) base_vars <- c("Age", base_vars)

  cat_vars <- c()
  if ("gender" %in% colnames(d))          cat_vars <- c(cat_vars, "gender")
  if ("pathologic_stage" %in% colnames(d)) cat_vars <- c(cat_vars, "pathologic_stage")
  if ("histological_type" %in% colnames(d)) cat_vars <- c(cat_vars, "histological_type")

  # keep categorical vars with >=2 levels
  cat_ok <- cat_vars[sapply(cat_vars, function(v) length(unique(na.omit(d[[v]]))) > 1)]

  # set histology reference if provided and present
  if ("histological_type" %in% cat_ok && !is.null(hist_ref_tbl)) {
    ref <- hist_ref_tbl %>% filter(acronym == cancer) %>% pull(hist_ref)
    if (length(ref) == 1 && ref %in% levels(d$histological_type)) {
      d$histological_type <- relevel(factor(d$histological_type), ref = ref)
    } else {
      d$histological_type <- factor(d$histological_type)
    }
  }

  vars_use <- c(base_vars, cat_ok)
  formula_str <- paste("Surv(days_to_censor, status) ~", paste(vars_use, collapse = " + "))
  fit <- coxph(as.formula(formula_str), data = d)

  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(acronym = cancer) %>%
    relocate(acronym)
}

# -----------------------------
# 3) Run across acronyms
# -----------------------------
cox_all <- df %>%
  group_by(acronym) %>%
  group_split() %>%
  map_df(~ run_cox_one_acronym(.x, hist_ref_tbl = hist_ref))

# Save full table
write_tsv(cox_all, OUT_ALL_TSV)

# HA-only summary (one row per acronym, if present)
cox_ha <- cox_all %>%
  filter(term == "HA.StrataHigh") %>%
  transmute(
    acronym,
    term,
    HR = estimate,
    conf.low,
    conf.high,
    p.value
  ) %>%
  arrange(p.value)

write_tsv(cox_ha, OUT_HA_TSV)

message("Saved:")
message("  ", OUT_ALL_TSV)
message("  ", OUT_HA_TSV)
