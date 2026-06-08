#!/bin/bash
set -euo pipefail

BASE_DIR="/home/sachleen.fnu/Projects/04282026_RNA_SEQ_Psomagen_Sachini"
ALIGNED="${BASE_DIR}/aligned"
GTF="${BASE_DIR}/genome/annotation.gtf"
COUNTS="${BASE_DIR}/counts"
THREADS=32
STRANDEDNESS=2   # 0=unstranded  1=forward  2=reverse (Zymo default)

source /home/sachleen.fnu/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq

mkdir -p "${COUNTS}"

# --- Strandedness verification (runs on 12M_a only, takes ~1 min) ---
echo "=== Strandedness verification on 12M_a ==="
TEST_BAM="${ALIGNED}/12M_a/Aligned.sortedByCoord.out.bam"

if [ ! -f "${TEST_BAM}" ]; then
    echo "ERROR: Test BAM not found: ${TEST_BAM}" >&2
    echo "Run 04_star_align.sh first." >&2
    exit 1
fi

for S in 0 1 2; do
    featureCounts -p -T 4 -s ${S} \
        -a "${GTF}" -t exon -g gene_id \
        -o /tmp/strand_test_s${S}.txt \
        "${TEST_BAM}" 2>/tmp/strand_test_s${S}.log
    ASSIGNED=$(grep "^Assigned" /tmp/strand_test_s${S}.txt.summary | awk '{print $2}')
    TOTAL=$(grep -v "^Status" /tmp/strand_test_s${S}.txt.summary | awk '{sum+=$2} END {print sum}')
    PCT=$(awk "BEGIN {printf \"%.1f\", ${ASSIGNED}/${TOTAL}*100}")
    echo "  -s ${S}: ${ASSIGNED} / ${TOTAL} reads assigned (${PCT}%)"
done

echo ""
echo "STRANDEDNESS in this script is set to: ${STRANDEDNESS}"
if [ -t 0 ]; then
    read -p "Press Enter to proceed, or Ctrl+C to abort and change STRANDEDNESS at top of script: "
else
    echo "Non-interactive mode: proceeding with STRANDEDNESS=${STRANDEDNESS}."
fi

# --- Collect all BAMs in sample order ---
mapfile -t BAMS < <(find "${ALIGNED}" -name "Aligned.sortedByCoord.out.bam" | sort)
if [ ${#BAMS[@]} -eq 0 ]; then
    echo "ERROR: No BAM files found under ${ALIGNED}" >&2
    exit 1
fi

echo "Running featureCounts on all 24 samples..."
featureCounts \
    -p \
    -T ${THREADS} \
    -s ${STRANDEDNESS} \
    -a "${GTF}" \
    -t exon \
    -g gene_id \
    -o "${COUNTS}/all_counts.txt" \
    "${BAMS[@]}"

echo "featureCounts complete."
echo "  Count matrix: ${COUNTS}/all_counts.txt"
echo "  Summary:      ${COUNTS}/all_counts.txt.summary"
