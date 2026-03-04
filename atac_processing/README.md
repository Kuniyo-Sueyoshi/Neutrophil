# ATAC-seq Processing Pipeline

## Overview

ATAC-seq data were processed against the human reference genome (GRCh38) using nf-core/atacseq.

The workflow includes:

- Adapter trimming
- Read alignment
- Peak calling
- Quality control (MultiQC)

All file paths are generalized and must be edited prior to execution.

---

## Directory Structure

atac_processing/
├── README.md
├── software_list.txt
├── data/
│   └── samplesheet.csv
├── script/
│   └── 01_nfcore_atacseq.sh
└── output/

---

## Software

Detailed software information is provided in `software_list.txt`.





