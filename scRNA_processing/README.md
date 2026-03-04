# scRNA Processing Pipeline

## Overview

This repository contains reproducible scripts for single-cell RNA-seq analyses performed in this study.

The workflow includes:

- FASTQ files processing (Cell Ranger)
- Quality control (QC) filtering  
- Ambient RNA evaluation (SoupX)  
- Doublet detection (scDblFinder)  
- Cross-species separation (human/mouse)  
- SCTransform-based integration (glmGamPoi)  
- Harmony batch correction  
- Differential expression and enrichment analyses  
- Human–mouse joint embedding  
- Ligand–receptor interaction analysis (NicheNet)

Full R session information (package versions and environment details) is provided in `sessionInfo.txt`.
All file paths are generalized and should be edited before execution.

---

## Directory structure

scRNA_processing/
├── README.md
├── sessionInfo.txt
├── R/
│   ├── fn.cluster.SCT.R
│   ├── fn.cluster.log.R
│   └── fn.QCmetrics.R
├── scripts/
│   ├── 00_cellranger_count.sh
│   ├── 01_Case01_preprocess.R
│   ├── 02_PDX_split_hg_mm.R
│   ├── 03_PDX_sample_integration_hg.R
│   ├── 04_PDX_sample_integration_mm.R
│   ├── 05_PDX_joint_embedding_hg_mm.R
│   └── 06_NicheNet_ligand_activity.R
├── data/
└── outp

---

## Software

Session Information is provided in `sessionInfo.txt` and `software_list.txt`.