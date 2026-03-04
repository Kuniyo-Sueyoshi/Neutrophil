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

```
wgs_processing/
├── README.md
├── software_versions.txt
├── data/
├── script/
│   ├── 01_gatk_bwa_hg38.sh
│   ├── 02_cnv_facets.sh
│   └── 03_mutect2.sh
└── output/
```

---


## Software

Software versions are listed in `software_versions.txt`.


