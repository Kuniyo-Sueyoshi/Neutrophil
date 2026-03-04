#!/usr/bin/env Rscript

############################################################
# 03_FAMD_and_Corr.HA.Neutrophil.R
#
# Purpose:
#   Explore relationships between HA signature and neutrophil fraction
#   in the MESOMICS cohort using:
#     1) FAMD (Factor Analysis for Mixed Data)
#     2) Correlation tests and scatter plot (HA vs Neutrophils)
#
# Inputs:
#   - "../output/MESOMICS_survival_table.rds"
#
# Outputs:
#   - Plots are created in the current R session (not written by default).
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(FactoMineR)
  library(factoextra)
  library(ggrepel)
})

# -----------------------------
# Paths (EDIT BEFORE RUN)
# -----------------------------
BASE_DIR <- "/PATH/TO/external_validation/MESOMICS"

IN_RDS <- file.path(BASE_DIR, "output", "MESOMICS_survival_table.rds")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1) Load annotation/features
# -----------------------------
dt.cc <- readRDS(IN_RDS)
dt.cc <- dt.cc %>%
  mutate(Type = Type %>% gsub("PM.B", "PM.BorS", .) %>% gsub("PM.S", "PM.BorS", .))
  
# -----------------------------
# 2) Prepare covariates for FAMD
# -----------------------------
# - Asbestos: combine professional + extra-professional exposure
# - Neutrophils: log-transform as in your original analysis
# - Type: collapse PM.B and PM.S -> PM.BorS
# FAMD input table
df.FAMD <- dt.cc %>%
  rename(HA.signature = HA) %>%
  rename(`%Neutrophil` = Neutrophils.log) %>%
  select(Age, Sex, Type, HA.signature, Asbestos, IHC.BAP1, `%Neutrophil`) %>%
  na.omit()

# -----------------------------
# 3) FAMD
# -----------------------------
res <- FAMD(df.FAMD, ncp = 5, graph = FALSE)

# variance explained (Dim1/Dim2)
variance1 <- round(res$eig[1, "percentage of variance"], 1)
variance2 <- round(res$eig[2, "percentage of variance"], 1)

# variable coordinates
var_type <- c(
  Age = "continuous",
  HA.signature = "continuous",
  Sex = "categorical",
  Type = "categorical",
  Asbestos = "categorical",
  IHC.BAP1 = "categorical",
  `%Neutrophil` = "continuous"
)

df.var <- res$var$coord %>%
  as.data.frame() %>%
  rownames_to_column("var") %>%
  select(var, Dim.1, Dim.2) %>%
  mutate(
    var.type = var_type[var],
    var = var %>% gsub("^Type$", "Hist.Type", .) %>% gsub("^IHC\\.BAP1$", "BAP1.status", .),
    oi = ifelse(var == "HA.signature", "HA.signature", "other"),
    var.type = factor(var.type, levels = c("continuous", "categorical"))
  )

p_famd <- ggplot(df.var, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(shape = var.type), size = 4) +
  ggrepel::geom_text_repel(aes(label = var), box.padding = 0.4) +
  theme_classic() +
  xlab(paste0("Dim1 (", variance1, "%)")) +
  ylab(paste0("Dim2 (", variance2, "%)"))

p_famd

# -----------------------------
# 4) Correlation: HA vs Neutrophils (overall and by histology)
# -----------------------------
cor_all <- cor.test(dt.cc$HA, dt.cc$Neutrophils.quanTIseq)
cor_E <- cor.test(
  dt.cc %>% filter(Type == "PM.E") %>% pull(HA),
  dt.cc %>% filter(Type == "PM.E") %>% pull(Neutrophils.quanTIseq)
)
cor_BorS <- cor.test(
  dt.cc %>% filter(Type == "PM.BorS") %>% pull(HA),
  dt.cc %>% filter(Type == "PM.BorS") %>% pull(Neutrophils.quanTIseq)
)

print(cor_all)
print(cor_E)
print(cor_BorS)

# -----------------------------
# 5) Scatter plot: HA vs %Neutrophil
# -----------------------------
p_scatter <- ggplot(dt.cc, aes(x = Neutrophils.quanTIseq, y = HA, color = Type)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  xlab("%Neutrophil (quanTIseq)") +
  ylab("HA signature")

p_scatter
