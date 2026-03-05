# External Validation and Survival Analysis Pipeline

## Overview

This repository contains analysis scripts used for survival modeling and validation of a histamine-response transcriptional signature in malignant pleural mesothelioma and pan-cancer cohorts.

The workflow includes two main components:

* **MESOMICS cohort analyses**
  Signature scoring, correlation analyses, survival modeling, and predictive model construction.

* **TCGA pan-cancer analyses**
  Immune deconvolution, external validation of the MESOMICS-trained model, and multivariable Cox analyses across cancer types.

The scripts are organized by cohort to separate **model development (MESOMICS)** from **external validation (TCGA)**.

---

# Directory Structure

```


# External Validation and Survival Analysis Pipeline

## Overview

This repository contains analysis scripts used for survival modeling and validation of a histamine-response transcriptional signature in malignant pleural mesothelioma and pan-cancer cohorts.

The workflow consists of two main components:

* **MESOMICS cohort analyses**
  Signature scoring, correlation analyses, survival modeling, and predictive model construction.

* **TCGA pan-cancer analyses**
  Immune deconvolution, external validation of the MESOMICS-trained model, and multivariable Cox analyses across cancer types.

The scripts are organized by cohort to clearly separate **model development (MESOMICS)** from **external validation (TCGA)**.

---

# Directory Structure

```
external_validation/
├── MESOMICS/
│   ├── data/
│   ├── output/
│   └── script/
│       ├── 01_Correlation_with_Case01.R
│       ├── 02_surv_HAsignature.R
│       ├── 03_FAMD_and_Corr.HA.Neutrophil.R
│       └── 04_surv.predict_model.construct.R
│
└── TCGA/
    ├── data/
    │   ├── lookup_hist_collapse.tsv
    │   ├── lookup_hist_reference.tsv
    │   └── lookup_stage_collapse.tsv
    │
    ├── output/
    │
    └── script/
        ├── 01_quantiseq_pancancer.R
        ├── 02_validate_lasso_cox_TCGA.R
        ├── 03_prepare_stage_histology.R
        └── 04_pancan_multivar.R
```

---

