#!/bin/bash
set -euo pipefail

BASE_DIR="/home/sachleen.fnu/Projects/04282026_RNA_SEQ_Psomagen_Sachini"
TRIMMED="${BASE_DIR}/trimmed"
STAR_INDEX="${BASE_DIR}/genome/star_index"
ALIGNED="${BASE_DIR}/aligned"
THREADS=32

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

mkdir -p "${ALIGNED}"

for SAMPLE in "${SAMPLES[@]}"; do
    echo "[$(date '+%H:%M:%S')] Aligning ${SAMPLE}..."
    mkdir -p "${ALIGNED}/${SAMPLE}"

    STAR \
        --runThreadN ${THREADS} \
        --genomeDir "${STAR_INDEX}" \
        --readFilesIn \
            "${TRIMMED}/${SAMPLE}_1_val_1.fq.gz" \
            "${TRIMMED}/${SAMPLE}_2_val_2.fq.gz" \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS NM \
        --outBAMsortingThreadN 16 \
        --outFileNamePrefix "${ALIGNED}/${SAMPLE}/"

    samtools index "${ALIGNED}/${SAMPLE}/Aligned.sortedByCoord.out.bam"

    echo "[$(date '+%H:%M:%S')] Done: ${SAMPLE}"
done

echo "STAR alignment complete. BAMs in: ${ALIGNED}"
