#!/bin/bash
set -euo pipefail

BASE_DIR="/home/sachleen.fnu/Projects/04282026_RNA_SEQ_Psomagen_Sachini"
RAW_DATA="${BASE_DIR}/raw_data"
TRIMMED="${BASE_DIR}/trimmed"
THREADS=16

SAMPLES=(
    12M_a  12M_b  12M_c  12M_d
    12R1_a 12R1_b 12R1_c 12R1_d
    12R25_a 12R25_b 12R25_c 12R25_d
    24M_a  24M_b  24M_c  24M_d
    24R1_a 24R1_b 24R1_c 24R1_d
    24R25_a 24R25_b 24R25_c 24R25_d
)

source /home/sachleen.fnu/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq

mkdir -p "${TRIMMED}"

for SAMPLE in "${SAMPLES[@]}"; do
    echo "[$(date '+%H:%M:%S')] Trimming ${SAMPLE}..."
    trim_galore \
        --paired \
        --cores ${THREADS} \
        -o "${TRIMMED}" \
        "${RAW_DATA}/${SAMPLE}_1.fastq.gz" \
        "${RAW_DATA}/${SAMPLE}_2.fastq.gz"
    echo "[$(date '+%H:%M:%S')] Done: ${SAMPLE}"
done

echo "Trim Galore complete. Trimmed reads in: ${TRIMMED}"
