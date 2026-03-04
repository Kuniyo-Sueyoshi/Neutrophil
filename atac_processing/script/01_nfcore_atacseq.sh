#!/usr/bin/env bash
set -euo pipefail

############################################################
# ATAC-seq processing using nf-core/atacseq
#
# Genome: GRCh38
# Read length: 150 bp
# Replicates: 2 (CellM)
############################################################

conda activate env_nf

NXF_OPTS='-Djava.io.tmpdir=/data/sueyoshi/tmp/' \
nextflow run nf-core/atacseq \
  --input ../data/samplesheet.csv \
  --outdir ../output/ \
  --genome GRCh38 \
  --read_length 150 \
  -profile docker
