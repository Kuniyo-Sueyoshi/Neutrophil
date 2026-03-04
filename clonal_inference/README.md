# Clonal Inference Pipeline

## Overview

Tumor clonal inference was performed by integrating whole-genome sequencing (WGS) and single-cell RNA sequencing (scRNA-seq) data from the same case.

Two complementary approaches were applied:

- PhyloWGS (bulk WGS–based clonal phylogeny)
- numbat (scRNA-seq allele-aware clonal inference)

All file paths are generalized and must be edited prior to execution.

---

## Directory Structure

```
clonal_inference/
├── README.md
├── sessionInfo.txt
├── software_list.txt
├── data/
├── script/
│   ├── 01_convert_facets_to_phylowgs.py
│   ├── 02_prepare_and_run_phylowgs.sh
│   ├── 03_pileup_phasing_numbat.sh
│   ├── 04_facets_to_seg_numbat.py
│   ├── 05_simplify_seg_numbat.py
│   └── 06_run_numbat.R
└── output/
```

---

## Workflow

### Bulk WGS (PhyloWGS)

- Somatic variant filtering
- FACETS CNV conversion
- PhyloWGS input generation
- MCMC-based clonal tree inference

### scRNA-seq (numbat)

- SNP pileup and phasing
- Optional WGS-derived CNV segmentation
- Clone inference using run_numbat()

---

## Software

Software versions are provided in `software_list.txt` and `sessionInfo.txt`.