#!/bin/bash
#
# prs_retry_failed.sh — Re-split failed PRS batches (005-010) into 250-ID chunks and resubmit
#
# The original batches at 500 IDs hit OOM (128GB) in FORMAT_SCOREFILES.
# This extracts the failed IDs and splits into smaller 250-ID batches.
#
# Usage: bash prs_retry_failed.sh [--dry-run]
#

set -euo pipefail

HPC_BASE="/rds/general/user/${HPC_USER:-YOUR_USERNAME}/home/GenomeAnalysis"
BATCH_DIR="${HPC_BASE}/results/prs/batches"
SCRIPTS_DIR="${HPC_BASE}/scripts"

DRY_RUN=false
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=true

# --- Collect failed batch IDs ---
FAILED_IDS="/tmp/failed_prs_ids.txt"
> "${FAILED_IDS}"

for batch_num in 005 006 007 008 009 010; do
    batch_file="${BATCH_DIR}/batch_${batch_num}.txt"
    if [[ -f "${batch_file}" ]]; then
        cat "${batch_file}" >> "${FAILED_IDS}"
        echo "  Collected: batch_${batch_num} ($(wc -l < "${batch_file}" | tr -d ' ') IDs)"
    else
        echo "  WARNING: ${batch_file} not found"
    fi
done

TOTAL_IDS=$(wc -l < "${FAILED_IDS}" | tr -d ' ')
echo ""
echo "Total failed IDs: ${TOTAL_IDS}"
echo "Splitting into 250-ID batches..."
echo ""

# --- Split into retry batches ---
RETRY_PREFIX="${BATCH_DIR}/batch_r"

# Clean any previous retry batches
rm -f "${BATCH_DIR}"/batch_r*.txt "${BATCH_DIR}"/batch_r*.done \
      "${BATCH_DIR}"/batch_r*.jobid "${BATCH_DIR}"/batch_r*.out \
      "${BATCH_DIR}"/batch_r*.err "${BATCH_DIR}"/batch_r*_job.sh

split -l 250 -d -a 3 "${FAILED_IDS}" "${RETRY_PREFIX}"

# Rename to .txt
for f in "${BATCH_DIR}"/batch_r[0-9][0-9][0-9]; do
    mv "$f" "${f}.txt"
done

N_BATCHES=$(ls "${BATCH_DIR}"/batch_r*.txt 2>/dev/null | wc -l | tr -d ' ')
echo "Created ${N_BATCHES} retry batches"
echo ""

# --- Submit each retry batch ---
for batch_file in "${BATCH_DIR}"/batch_r*.txt; do
    BATCH_NAME=$(basename "$batch_file" .txt)
    N_IN_BATCH=$(wc -l < "$batch_file" | tr -d ' ')

    echo "Batch: ${BATCH_NAME} (${N_IN_BATCH} IDs)"

    if [[ "$DRY_RUN" == true ]]; then
        echo "  [DRY RUN] Would submit with 128GB"
        continue
    fi

    JOB_SCRIPT="${BATCH_DIR}/${BATCH_NAME}_job.sh"
    cat > "$JOB_SCRIPT" << PBSEOF
#!/bin/bash
#PBS -N prs_${BATCH_NAME}
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=8:mem=128gb
#PBS -o ${BATCH_DIR}/${BATCH_NAME}.out
#PBS -e ${BATCH_DIR}/${BATCH_NAME}.err
#PBS -j n

set -euo pipefail

# cd to per-batch output dir to avoid .nextflow/ lock contention
BATCH_OUTDIR="${HPC_BASE}/results/prs/${BATCH_NAME}"
mkdir -p "\${BATCH_OUTDIR}"
cd "\${BATCH_OUTDIR}"

bash ${SCRIPTS_DIR}/run_pgscalc.sh --file ${batch_file} --batch-name ${BATCH_NAME}

touch "${BATCH_DIR}/${BATCH_NAME}.done"
PBSEOF

    JOB_ID=$(/opt/pbs/bin/qsub "$JOB_SCRIPT")
    echo "  Submitted: ${JOB_ID}"
    echo "${JOB_ID}" > "${BATCH_DIR}/${BATCH_NAME}.jobid"
    sleep 2
done

echo ""
echo "All retry batches submitted."
echo "Monitor: bash prs_batch.sh --status"
