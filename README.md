# analysis

## Overview

This directory contains the analysis scripts used in the study:

**“Cooperative Response Loops Enhanced by Histamine Signaling during Tumor–Neutrophil Crosstalk in Pleural Mesothelioma.”**

The repository provides preprocessing and analysis pipelines for sequencing-based experiments used in the manuscript, including bulk RNA-seq, scRNA-seq, ATAC-seq, WGS, and external validation analyses.

## Directory structure

```
analysis/
├── atac_processing/        ATAC-seq preprocessing and analysis
├── bulkRNA_processing/     Bulk RNA-seq processing and differential expression
├── clonal_inference/       Tumor clonal inference (FACETS / PhyloWGS / Numbat)
├── cnv_integration/        Integration of CNV segments across datasets
├── scRNA_processing/       Single-cell RNA-seq preprocessing and analysis
├── wgs_processing/         Whole-genome sequencing preprocessing
└── external_validation/    Independent validation analyses (MESOMICS and TCGA)
```

## Input dataset

Input datasets (e.g., sequencing reads and large expression matrices) are not included in this repository and should be obtained from the corresponding public repositories or project data sources, stated in the manuscript.

## Requirements

Analyses were performed using a combination of:

* **R** (Seurat, GSVA, survival, glmnet, tidyverse)
* **Python**
* **Shell pipelines**
* External bioinformatics tools (e.g., GATK, FACETS, PhyloWGS, Numbat)

Exact tool usage is documented within individual directories.

## Notes

Paths to external datasets are intentionally generalized in the scripts and should be adapted to the local environment before execution.
