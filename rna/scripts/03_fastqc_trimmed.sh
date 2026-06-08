#!/bin/bash
set -euo pipefail

BASE_DIR="/home/sachleen.fnu/Projects/04282026_RNA_SEQ_Psomagen_Sachini"
TRIMMED="${BASE_DIR}/trimmed"
QC_TRIMMED="${BASE_DIR}/qc/fastqc_trimmed"
THREADS=32

source /home/sachleen.fnu/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq

mkdir -p "${QC_TRIMMED}"

echo "Running FastQC on trimmed reads..."
fastqc \
    -t ${THREADS} \
    -o "${QC_TRIMMED}" \
    "${TRIMMED}"/*_val_1.fq.gz \
    "${TRIMMED}"/*_val_2.fq.gz

echo "FastQC (trimmed) complete. Reports in: ${QC_TRIMMED}"
