#!/bin/bash
set -euo pipefail

BASE_DIR="/home/sachleen.fnu/Projects/04282026_RNA_SEQ_Psomagen_Sachini"

mkdir -p \
    "${BASE_DIR}/raw_data" \
    "${BASE_DIR}/genome/star_index" \
    "${BASE_DIR}/qc/fastqc_raw" \
    "${BASE_DIR}/qc/fastqc_trimmed" \
    "${BASE_DIR}/qc/multiqc" \
    "${BASE_DIR}/trimmed" \
    "${BASE_DIR}/aligned" \
    "${BASE_DIR}/counts" \
    "${BASE_DIR}/scripts"

echo "Directory structure created under ${BASE_DIR}"
echo ""
echo "Before running the pipeline, place:"
echo "  - FASTQ files in:       ${BASE_DIR}/raw_data/  (named <sample>_1.fastq.gz and <sample>_2.fastq.gz)"
echo "  - STAR index in:        ${BASE_DIR}/genome/star_index/"
echo "  - GTF annotation at:    ${BASE_DIR}/genome/annotation.gtf"

