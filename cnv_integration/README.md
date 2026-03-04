# CNV Integration Pipeline

## Overview

Copy number variation (CNV) signals derived from bulk whole-genome sequencing (WGS) and single-cell RNA-seq (inferCNV) were integrated using a window-segmented genome representation.

The workflow includes:

- Genome segmentation into fixed windows
- Conversion of FACETS CNV segments to segmented windows
- Aggregation of CNV signals across segmented genome bins
- Correlation analysis between WGS-derived CNV and inferCNV signals
- Projection of inferCNV cell-group signals to single-cell resolution

File paths are generalized and must be edited prior to execution.

---

## Directory Structure

```
cnv_integration/
├── README.md
├── sessionInfo.txt
├── software_versions.txt
├── data/
│   └── chr_len.tsv
├── output/
│   └── chr_segmented.tsv
└── script/
    ├── 00_segment_template.py
    ├── 01_segment_cnv.py
    ├── 02_sc_cnv_correlation.py
    └── utils_make_SegChr.py
```

---

## Software

Software versions are documented in `software_versions.txt`.
