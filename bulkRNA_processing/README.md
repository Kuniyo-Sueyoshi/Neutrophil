# Bulk RNA-seq Processing Pipeline

## Overview

Bulk RNA-seq differential expression analyses were performed in R (v4.2.3) using edgeR and limma-voom.

The workflow includes:

- Transcript-level quantification (Salmon)
- Gene-level summarization (tximport)
- edgeR normalization
- limma-voom differential expression
- GSVA pathway analysis (neutrophils only)

File paths are generalized and must be edited prior to execution.

---

## Directory Structure

```
bulkRNA_processing/
├── README.md
├── sessionInfo.txt
├── software_list.txt
├── data/
├── output/
└── script/
    ├── 01_salmon_quant.sh
    ├── 02_tximport_summarize.R
    ├── 03_CellM_bulk_limma_voom.R
    └── 04_Neutro_bulk_limma_voom.R    

```

---

## Software

software versions are documented in `sessionInfo.txt` and `software_list.txt`.




