# WGS Processing Pipeline

## Overview

Whole-genome sequencing (WGS) analyses were performed against the human reference genome (hg38).

The workflow includes:

- Read alignment (BWA-MEM)
- BAM sorting and indexing (samtools)
- Duplicate marking and base quality score recalibration (GATK)
- Allele-specific copy number analysis (FACETS)
- Somatic variant calling (GATK Mutect2)

All file paths are generalized and must be edited prior to execution.

---

## Directory Structure

wgs_processing/
├── README.md
├── software_versions.txt
├── data/
├── script/
│   ├── 01_gatk_bwa_hg38.sh
│   ├── 02_cnv_facets.sh
│   └── 03_mutect2.sh
└── output/

---

## Pipeline Steps

### 1. Alignment and Preprocessing
`01_gatk_bwa_hg38.sh`

- BWA-MEM alignment
- Sorting and indexing
- MarkDuplicates
- Base Quality Score Recalibration (BQSR)

### 2. Copy Number Analysis
`02_cnv_facets.sh`

- Allele-specific CNV analysis using FACETS

### 3. Somatic Variant Calling
`03_mutect2.sh`

- GATK Mutect2
- FilterMutectCalls

---

## Software

Software versions are listed in `software_versions.txt`.


