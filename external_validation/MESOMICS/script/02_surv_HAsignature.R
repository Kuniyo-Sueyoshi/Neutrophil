#!/usr/bin/env Rscript

############################################################
# 02_surv_HAsignature.R
#
# Purpose:
#   Survival analysis in MESOMICS using a histamine-response signature
#   derived from CellM (DEG: HA vs DW).
#
# Inputs (generalized):
#   - MESOMICS annotation table: data/MESOMICS_ann_110pt.tsv
#   - MESOMICS gene TPM matrix (Symbol): data/gene_TPM_matrix_1pass_Symbol.tsv
#   - CellM DEG table (HA vs DW): /PATH/DEG_CellM_HA_vs_DW.tsv
#
# Analysis:
#   1) GSVA scoring (HA gene set; log2(TPM+1))
#   2) Kaplan–Meier (Q4 vs Q1-3)
#   3) Cox regression (univariate and multivariate)
#
# Outputs:
#   - output/MESOMICS_survival_table.rds
#   - (interactive) KM plot object: g3
#   - (interactive) Cox tables: univ_tbl, multi_tbl
# 
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(survminer)
  library(GSVA)
  library(broom)
})

# -----------------------------
# Paths (EDIT BEFORE RUN)
# -----------------------------
BASE_DIR <- "/PATH/TO/external_validation/MESOMICS"

IN_ANN   <- file.path(BASE_DIR, "data", "MESOMICS_ann_110pt.tsv")
IN_TPM   <- file.path(BASE_DIR, "data", "gene_TPM_matrix_1pass_Symbol.tsv")
IN_DEG   <- "/PATH/DEG_CellM_HA_vs_DW.tsv"
OUT_DIR <- file.path(BASE_DIR, "output")

# -----------------------------
# 1) Load inputs
# -----------------------------
ann <- read_tsv(IN_ANN, show_col_types = FALSE)

d.logTPM <- read_tsv(IN_TPM, show_col_types = FALSE) %>%
  column_to_rownames("Symbol") %>%
  { log2(. + 1) }

# Align samples
comm <- intersect(ann$ID, colnames(d.logTPM))
ann <- ann[match(comm, ann$ID), , drop = FALSE]
d.logTPM <- d.logTPM[, match(comm, colnames(d.logTPM)), drop = FALSE]
stopifnot(all(ann$ID == colnames(d.logTPM)))

# Histamine-response gene set (CellM HA-up)
gene.ha <- read_tsv(IN_DEG, show_col_types = FALSE) %>%
  filter(logFC > 0, adj.P.Val < 0.05) %>%
  pull(Gene)

# -----------------------------
# 2) GSVA score (single gene set)
# -----------------------------
gsva.es <- gsva(
  gsvaParam(
    exprData = as.matrix(d.logTPM),
    geneSets = list(HA = gene.ha),
    kcdf = "Gaussian"
  )
)

dt.c <- merge(
  t(gsva.es) %>% as.data.frame() %>% rownames_to_column("ID"),
  ann,
  by = "ID"
) %>%
  mutate(
    Survival.Time = Survival.Time / 12,                       # months -> years
    Survival.Censor.digit = ifelse(Survival.Censor == "dead", 1, 0)  # event=1
  )

# -----------------------------
# 3) KM plot (Q4 vs Q1-3) : HA
# -----------------------------
l.ha <- dt.c$HA > quantile(dt.c$HA, probs = 0.75, na.rm = TRUE)

dt.c1 <- dt.c %>%
  mutate(Strata = ifelse(l.ha, "Q4", "Q1-3"))

fit1 <- survfit(Surv(Survival.Time, Survival.Censor.digit) ~ Strata, data = dt.c1)

g3 <- ggsurvplot(
  fit1,
  data = dt.c1,
  pval = TRUE,
  risk.table = TRUE,
  break.x.by = 2,
  surv.median.line = "h",
  pval.coord = c(3.5, 0.4)
) +
  ggtitle("MESOMICS cohort") +
  xlab("Time (year)")


# ============================================================
# 5) Cox regression (univariate / multivariate)
# ============================================================

# Assemble model table
dt.cc <- dt.c %>%
  mutate(
    HA.Strata = ifelse(l.ha, "Q4", "Q1-3"),
    Type = Type %>% gsub("PM.B", "PM.BorS", .) %>% gsub("PM.S", "PM.BorS", .),
    Type = factor(Type, levels = c("PM.E", "PM.BorS"))
  ) %>%
  select(
    Age, Sex, Type, HA.Strata, HA, CXCL5, Neutrophils.quanTIseq, IHC.BAP1,
    Survival.Censor.digit, Survival.Time,
    Professional.Asbestos, ExtraProfessional.Asbestos
  ) %>%
  mutate(
    Asbestos = Professional.Asbestos,
    Asbestos = ifelse(ExtraProfessional.Asbestos == "Exposed", "Exposed", Asbestos),
    Neutrophils.log = log(Neutrophils.quanTIseq * 100 + 1)
  )

# -----------------------------
# Univariate Cox
# -----------------------------
vars_uni <- c("Age", "Sex", "Type", "Asbestos", "HA.Strata")

univ_fits <- lapply(vars_uni, function(v) {
  f <- as.formula(paste0("Surv(Survival.Time, Survival.Censor.digit) ~ ", v))
  coxph(f, data = dt.cc)
})
names(univ_fits) <- vars_uni

univ_tbl <- bind_rows(lapply(names(univ_fits), function(v) {
  broom::tidy(univ_fits[[v]], exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(Variable = v)
})) %>%
  mutate(
    ReferenceLevel = ifelse(Variable == "Age", "-", NA_character_),
    ComparedLevel  = ifelse(Variable == "Age", "-", term),
    `HR [95%CI]`   = sprintf("%.2f [%.2f–%.2f]", estimate, conf.low, conf.high)
  ) %>%
  select(Variable, ReferenceLevel, ComparedLevel, `HR [95%CI]`, p.value)

univ_tbl

# -----------------------------
# Multivariate Cox
# -----------------------------
# Keep this minimal; extend as needed (e.g., + Sex + Asbestos + IHC.BAP1)
cox_multi <- coxph(
  Surv(Survival.Time, Survival.Censor.digit) ~ Age + Type + HA.Strata,
  data = dt.cc
)

multi_tbl <- broom::tidy(cox_multi, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    Variable = case_when(
      term == "Age" ~ "Age",
      str_starts(term, "Type") ~ "Type",
      str_starts(term, "HA.Strata") ~ "HA.Strata",
      TRUE ~ term
    ),
    ReferenceLevel = case_when(
      Variable == "Age" ~ "-",
      Variable == "Type" ~ levels(dt.cc$Type)[1],
      Variable == "HA.Strata" ~ levels(factor(dt.cc$HA.Strata))[1],
      TRUE ~ "-"
    ),
    ComparedLevel = case_when(
      Variable == "Age" ~ "-",
      TRUE ~ term
    ),
    `HR [95%CI]` = sprintf("%.2f [%.2f–%.2f]", estimate, conf.low, conf.high)
  ) %>%
  select(Variable, ReferenceLevel, ComparedLevel, `HR [95%CI]`, p.value)

multi_tbl


# --------------------------------------------------
# Save analysis table for downstream analyses
# --------------------------------------------------
dir.create(OUT_DIR, showWarnings = FALSE)
saveRDS(dt.cc, file = file.path(OUT_DIR, "MESOMICS_survival_table.rds"))