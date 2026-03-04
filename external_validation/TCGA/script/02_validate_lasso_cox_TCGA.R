#!/usr/bin/env Rscript

############################################################
# 02_validate_lasso_cox_TCGA.R
#
# Purpose:
#   Validate the MESOMICS-trained Cox LASSO model in TCGA-MESO.
#   - Compute HA signature (GSVA) and extract CXCL5 expression
#   - Merge with TCGA clinical + BAP1 status + quanTIseq neutrophils
#   - Build design matrix with the same covariates as training
#   - Predict risk score and evaluate:
#       * C-index
#       * KM (High vs Low; median split)
#
# Inputs (generalized):
#   - MESOMICS model (trained): ../MESOMICS/output/lasso_cox_model.rds
#   - TCGA expression (gene x sample or sample x gene):
#   - DEG table for HA signature (CellM; HA vs DW): See "bulkRNA_processing"
#   - TCGA clinical annotation
#   - TCGA BAP1 mutation/CNV status (from cBioPortal)
#   - quanTIseq output (neutrophils): (prepare separately)
#
# Outputs:
#   - output/validation_table_TCGA_MESO.tsv
#   - output/KM_TCGA_validation.pdf
#
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(survminer)
  library(glmnet)
  library(GSVA)
})

set.seed(123)

# -----------------------------
# Paths (EDIT BEFORE RUN)
# -----------------------------
BASE_DIR_TCGA <- "/PATH/TO/external_validation/TCGA"
BASE_DIR_MESO <- "/PATH/TO/external_validation/MESOMICS"

IN_MODEL_RDS <- file.path(BASE_DIR_MESO, "output", "lasso_cox_model.rds")

IN_EXPR_TSV  <- file.path(BASE_DIR_TCGA, "data", "TCGA_RNAseqV2_Meso_RSEM_TPM.tsv")
IN_CLIN_TSV  <- file.path(BASE_DIR_TCGA, "data", "clinical_annotation_TCGA_MESO.tsv")
IN_BAP1_TSV  <- file.path(BASE_DIR_TCGA, "data", "TCGA.mutation_BAP1_status_.tsv")

IN_DEG_TSV   <- "/PATH/DEG_CellM_HA_vs_DW.tsv"

# quanTIseq results 
IN_QUANTISEQ_TSV <- file.path(BASE_DIR_TCGA, "output", "Immune10_deconvolv_TCGA.PanCan.annotated.tsv")


# Output
OUT_VAL_TSV <- file.path(BASE_DIR_TCGA, "output", "validation_table_TCGA_MESO.tsv")
OUT_KM_PDF <- file.path(BASE_DIR_TCGA, "output", "KM_TCGA_validation.pdf")

# -----------------------------
# 0) Load trained model
# -----------------------------
fit.lasso <- readRDS(IN_MODEL_RDS)

# -----------------------------
# 1) Load TCGA expression and log-transform
# -----------------------------
# Expected: rows = genes (Symbol), columns = samples (bcr_patient_barcode)
# If your table is the opposite, transpose after reading.
d.TPM.TCGA <- read_tsv(IN_EXPR_TSV, show_col_types = FALSE)
d.logTPM.TCGA <- log2(d.TPM.TCGA + 1)

# -----------------------------
# 2) Build HA gene set from CellM DEG
# -----------------------------
deg_cellm <- read_tsv(IN_DEG_TSV, show_col_types = FALSE)
gene.ha <- deg_cellm %>%
  filter(logFC > 0, adj.P.Val < 0.05) %>%
  pull(Gene) %>%
  unique()

# -----------------------------
# 3) HA signature by GSVA + CXCL5 extraction
# -----------------------------
gsva.es <- gsva(
  gsvaParam(
    exprData = as.matrix(d.logTPM.TCGA),
    geneSets = list(HA = gene.ha),
    kcdf = "Gaussian"
  )
)

# Ensure CXCL5 exists
stopifnot("CXCL5" %in% rownames(d.logTPM.TCGA))

dt.ha.cxcl5 <- rbind(
  CXCL5 = d.logTPM.TCGA["CXCL5", ],
  HA    = gsva.es["HA", ]
) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("bcr_patient_barcode")

# -----------------------------
# 4) Load annotations (clinical, BAP1, quanTIseq)
# -----------------------------
meta.data <- read_tsv(IN_CLIN_TSV, show_col_types = FALSE)
mut.data  <- read_tsv(IN_BAP1_TSV, show_col_types = FALSE)

quantISeq.data.pancan <- read_tsv(IN_QUANTISEQ_TSV, show_col_types = FALSE)
quantISeq.data <- quantISeq.data.pancan %>%
  filter(acronym == "MESO") %>%
  select(Sample, Neutrophils)

# -----------------------------
# 5) Merge into validation table (dt.val)
# -----------------------------
anno <- merge(
  meta.data %>%
    select(
      bcr_patient_barcode,
      age_at_initial_pathologic_diagnosis,
      gender,
      histological_type,
      days_to_death,
      vital_status,
      status,
      days_to_censor
    ),
  mut.data %>%
    select(Patient_ID, BAP1_Altered, BAP1_Detail, Significance),
  by.x = "bcr_patient_barcode", by.y = "Patient_ID",
  all.x = TRUE
) %>%
  tidyr::replace_na(replace = list(BAP1_Altered = 0, Significance = "wt")) %>%
  mutate(IHC.BAP1 = ifelse(BAP1_Altered == 1, "Loss", "Retain")) %>%
  merge(
    .,
    quantISeq.data %>% select(Sample, Neutrophils),
    by.x = "bcr_patient_barcode", by.y = "Sample",
    all.x = TRUE
  ) %>%
  merge(
    .,
    dt.ha.cxcl5,
    by = "bcr_patient_barcode",
    all.x = TRUE
  ) %>%
  mutate(
    Survival.Time = as.numeric(days_to_censor) / 365,
    Sex = ifelse(gender == "MALE", "M", "F"),
    Type = histological_type %>%
      gsub("MM.nos", "PM.E", .) %>%
      gsub("MM", "PM", .),
    # if TCGA has PM.S / PM.B, collapse to PM.BorS
    Type = Type %>% gsub("PM.B", "PM.BorS", .) %>% gsub("PM.S", "PM.BorS", .),
    Type = factor(Type, levels = c("PM.E", "PM.BorS")),
    Age = as.numeric(age_at_initial_pathologic_diagnosis),
    Survival.Censor.digit = as.numeric(status),
    Neutrophils.log = log(Neutrophils * 100 + 1)
  ) %>%
  dplyr::rename(HA = HA)

dt.val <- anno %>%
  filter(!is.na(Survival.Time), !is.na(Survival.Censor.digit)) %>%
  mutate(Survival.Time = ifelse(Survival.Time <= 0, 1e-6, Survival.Time))

# -----------------------------
# 6) Build design matrix for validation
# -----------------------------
x_val <- model.matrix(
  ~ Age + Sex + Type + HA + CXCL5 + Neutrophils.log + IHC.BAP1,
  data = dt.val
)[, -1, drop = FALSE]

# -----------------------------
# 7) Predict risk score + evaluate
# -----------------------------
risk_score_full <- predict(fit.lasso, newx = x_val, s = "lambda.min", type = "link")
dt.val$risk_score_full <- as.numeric(risk_score_full)

dt.val$risk_group_full <- ifelse(
  dt.val$risk_score_full > median(dt.val$risk_score_full, na.rm = TRUE),
  "High", "Low"
)
dt.val$risk_group_full <- factor(dt.val$risk_group_full, levels = c("Low", "High"))

fit.km <- survfit(
  Surv(Survival.Time, Survival.Censor.digit) ~ risk_group_full,
  data = dt.val
)

c_index_val <- survival::concordance(
  Surv(Survival.Time, Survival.Censor.digit) ~ risk_score_full,
  data = dt.val,
  reverse = TRUE
)$concordance

message("C-index (TCGA): ", round(c_index_val, 3))

g <- ggsurvplot(
  fit.km,
  data = dt.val,
  pval = TRUE,
  risk.table = TRUE,
  surv.median.line = "h",
  pval.coord = c(3, 0.5),
  title = "Validation cohort (TCGA)",
  xlab = "Time (year)",
  legend.title = "Risk group",
  legend.labs = c("Low", "High")
)

ggsave(OUT_KM_PDF, g$plot, width = 7, height = 6)
write_tsv(dt.val, OUT_VAL_TSV)