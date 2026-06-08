#!/bin/bash
set -euo pipefail

BASE_DIR="/home/sachleen.fnu/Projects/04282026_RNA_SEQ_Psomagen_Sachini"
RAW_DATA="${BASE_DIR}/raw_data"
QC_RAW="${BASE_DIR}/qc/fastqc_raw"
THREADS=32

source /home/sachleen.fnu/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq

mkdir -p "${QC_RAW}"

echo "Running FastQC on raw reads..."
fastqc \
    -t ${THREADS} \
    -o "${QC_RAW}" \
    "${RAW_DATA}"/*_1.fastq.gz \
    "${RAW_DATA}"/*_2.fastq.gz

echo "FastQC (raw) complete. Reports in: ${QC_RAW}"
