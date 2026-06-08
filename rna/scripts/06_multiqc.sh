#!/bin/bash
set -euo pipefail

BASE_DIR="/home/sachleen.fnu/Projects/04282026_RNA_SEQ_Psomagen_Sachini"
MULTIQC_OUT="${BASE_DIR}/qc/multiqc"

source /home/sachleen.fnu/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq

mkdir -p "${MULTIQC_OUT}"

echo "Running MultiQC..."
multiqc \
    "${BASE_DIR}/qc/fastqc_raw/" \
    "${BASE_DIR}/qc/fastqc_trimmed/" \
    "${BASE_DIR}/trimmed/" \
    "${BASE_DIR}/aligned/" \
    "${BASE_DIR}/counts/" \
    --outdir "${MULTIQC_OUT}" \
    --force

echo "MultiQC complete. Report: ${MULTIQC_OUT}/multiqc_report.html"
