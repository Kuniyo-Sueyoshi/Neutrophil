#!/usr/bin/env Rscript

############################################################
# 04_surv.predict_model.construct.R
#
# Purpose:
#   Construct a Cox LASSO survival model (glmnet) in MESOMICS
#   and evaluate KM by median risk score.
#
# Input:
#   - output/MESOMICS_survival_table.rds
#
# Output:
#   - output/lasso_cox_model.rds
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(glmnet)
  library(survminer)
})

set.seed(123)

# -----------------------------
# Paths
# -----------------------------
BASE_DIR <- "/PATH/TO/external_validation/MESOMICS"

IN_RDS <- file.path(BASE_DIR, "output", "MESOMICS_survival_table.rds")
OUT_RDS <- file.path(BASE_DIR, "output", "lasso_cox_model.rds")

# -----------------------------
# Load data
# -----------------------------
dt.cc <- readRDS(IN_RDS)

# -----------------------------
# Prepare modeling table
# -----------------------------
dt.glm <- dt.cc %>%
  select(Survival.Time, Survival.Censor.digit, Age, Sex, Type, HA, CXCL5, Neutrophils.quanTIseq, IHC.BAP1) %>%
  drop_na()
dt.glm$Survival.Time[dt.glm$Survival.Time <= 0] <- 1e-6
dt.glm <- dt.glm %>%
  mutate(
    Neutrophils.log = log(Neutrophils.quanTIseq * 100 + 1)
  )

# -----------------------------
# Design matrix
# -----------------------------
x <- model.matrix(
  Surv(Survival.Time, Survival.Censor.digit) ~
    Age + Sex + Type + HA + CXCL5 + Neutrophils.log + IHC.BAP1,
  data = dt.glm
)[, -1]

y <- with(dt.glm, Surv(Survival.Time, Survival.Censor.digit))

# -----------------------------
# LASSO Cox model
# -----------------------------
fit.lasso <- cv.glmnet(
  x = x,
  y = y,
  family = "cox",
  nfolds = 5,
  type.measure = "C",
  alpha = 1
)

saveRDS(fit.lasso, OUT_RDS)

# -----------------------------
# Risk score
# -----------------------------
dt.glm$risk_score <- as.numeric(
  predict(fit.lasso, newx = x, s = "lambda.min", type = "link")
)

# -----------------------------
# C-index
# -----------------------------
c_index_val <- survival::concordance(
  Surv(Survival.Time, Survival.Censor.digit) ~ risk_score,
  data = dt.glm,
  reverse = TRUE
)$concordance

message("C-index (apparent): ", round(c_index_val, 3))

# -----------------------------
# KM by median risk score
# -----------------------------
dt.glm$risk_group <- ifelse(
  dt.glm$risk_score > median(dt.glm$risk_score),
  "High", "Low"
)
dt.glm$risk_group <- factor(dt.glm$risk_group, levels = c("Low","High"))

fit.km <- survfit(
  Surv(Survival.Time, Survival.Censor.digit) ~ risk_group,
  data = dt.glm
)

g <- ggsurvplot(
  fit.km,
  data = dt.glm,
  pval = TRUE,
  risk.table = TRUE,
  surv.median.line = "h",
  title = "Training cohort (MESOMICS)",
  xlab = "Time (year)"
)

g
