#!/bin/bash
#
# prs_batch.sh — Orchestrate multi-batch pgsc_calc runs
#
# Splits a large list of PGS IDs into batches and submits each as a
# separate pgsc_calc run on the HPC. Concatenates results after all
# batches complete.
#
# Usage:
#   # Step 1: Download all PGS IDs (run on login node — needs internet)
#   python3 download_pgs_ids.py --output pgs_ids.txt --metadata --metadata-output pgs_metadata.json
#
#   # Step 2: Split into batches and submit
#   bash prs_batch.sh --ids pgs_ids.txt [--batch-size 500] [--dry-run]
#
#   # Step 3: Check status
#   bash prs_batch.sh --status
#
#   # Step 4: Concatenate results when all done
#   bash prs_batch.sh --concat
#

set -euo pipefail

# --- Configuration ---
HPC_BASE="/rds/general/user/${HPC_USER:-YOUR_USERNAME}/home/GenomeAnalysis"
RESULTS_DIR="${HPC_BASE}/results/prs"
BATCH_DIR="${RESULTS_DIR}/batches"
SCRIPTS_DIR="${HPC_BASE}/scripts"
BATCH_SIZE=500
DRY_RUN=false

# --- Parse arguments ---
IDS_FILE=""
ACTION=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --ids|-i)
            IDS_FILE="$2"
            ACTION="submit"
            shift 2
            ;;
        --batch-size|-b)
            BATCH_SIZE="$2"
            shift 2
            ;;
        --status|-s)
            ACTION="status"
            shift
            ;;
        --concat|-c)
            ACTION="concat"
            shift
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

if [[ -z "$ACTION" ]]; then
    echo "Usage:"
    echo "  $0 --ids pgs_ids.txt [--batch-size 500] [--dry-run]"
    echo "  $0 --status"
    echo "  $0 --concat"
    exit 1
fi

# --- Submit batches ---
if [[ "$ACTION" == "submit" ]]; then
    if [[ -z "$IDS_FILE" || ! -f "$IDS_FILE" ]]; then
        echo "ERROR: PGS ID file not found: $IDS_FILE" >&2
        exit 1
    fi

    TOTAL_IDS=$(wc -l < "$IDS_FILE" | tr -d ' ')
    N_BATCHES=$(( (TOTAL_IDS + BATCH_SIZE - 1) / BATCH_SIZE ))

    echo "=== PRS Batch Submission ==="
    echo "Total PGS IDs: ${TOTAL_IDS}"
    echo "Batch size: ${BATCH_SIZE}"
    echo "Number of batches: ${N_BATCHES}"
    echo ""

    # Clean batch directory from previous runs
    rm -rf "${BATCH_DIR}"
    mkdir -p "${BATCH_DIR}"

    # Split IDs into batch files (creates batch_000, batch_001, etc.)
    split -l "${BATCH_SIZE}" -d -a 3 "$IDS_FILE" "${BATCH_DIR}/batch_"

    # Rename split output to .txt extension
    # split creates files without extension — only rename those (no dots in name)
    for f in "${BATCH_DIR}"/batch_[0-9][0-9][0-9]; do
        mv "$f" "${f}.txt"
    done

    # Submit each batch
    BATCH_NUM=0
    for batch_file in "${BATCH_DIR}"/batch_*.txt; do
        BATCH_NUM=$((BATCH_NUM + 1))
        BATCH_NAME=$(basename "$batch_file" .txt)
        N_IN_BATCH=$(wc -l < "$batch_file" | tr -d ' ')

        echo "Batch ${BATCH_NUM}/${N_BATCHES}: ${BATCH_NAME} (${N_IN_BATCH} scores)"

        if [[ "$DRY_RUN" == true ]]; then
            echo "  [DRY RUN] Would submit: run_pgscalc.sh --file ${batch_file}"
        else
            # Create a PBS job script for this batch
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

# cd to per-batch output dir to avoid .nextflow/ lock contention in shared CWD
BATCH_OUTDIR="${HPC_BASE}/results/prs/${BATCH_NAME}"
mkdir -p "\${BATCH_OUTDIR}"
cd "\${BATCH_OUTDIR}"

# Run pgsc_calc for this batch (no --resume on fresh run)
bash ${SCRIPTS_DIR}/run_pgscalc.sh --file ${batch_file} --batch-name ${BATCH_NAME}

# Mark batch as complete (only reached if pgsc_calc succeeds)
touch "${BATCH_DIR}/${BATCH_NAME}.done"
PBSEOF

            # Submit the PBS job
            JOB_ID=$(/opt/pbs/bin/qsub "$JOB_SCRIPT")
            echo "  Submitted: ${JOB_ID}"
            echo "${JOB_ID}" > "${BATCH_DIR}/${BATCH_NAME}.jobid"

            # Brief pause between submissions
            sleep 2
        fi
    done

    echo ""
    echo "All batches submitted. Monitor with: $0 --status"
fi

# --- Check status ---
if [[ "$ACTION" == "status" ]]; then
    echo "=== PRS Batch Status ==="

    if [[ ! -d "$BATCH_DIR" ]]; then
        echo "No batches found. Run with --ids first."
        exit 0
    fi

    TOTAL=0
    DONE=0
    RUNNING=0
    FAILED=0

    for batch_file in "${BATCH_DIR}"/batch_*.txt; do
        BATCH_NAME=$(basename "$batch_file" .txt)
        TOTAL=$((TOTAL + 1))

        if [[ -f "${BATCH_DIR}/${BATCH_NAME}.done" ]]; then
            DONE=$((DONE + 1))
            echo "  [DONE]    ${BATCH_NAME}"
        elif [[ -f "${BATCH_DIR}/${BATCH_NAME}.jobid" ]]; then
            JOB_ID=$(cat "${BATCH_DIR}/${BATCH_NAME}.jobid")
            # Check if job is still running
            if /opt/pbs/bin/qstat "$JOB_ID" 2>/dev/null | grep -q "[RQ]"; then
                RUNNING=$((RUNNING + 1))
                echo "  [RUNNING] ${BATCH_NAME} (${JOB_ID})"
            else
                # Job finished but no .done marker — check for errors
                if [[ -f "${BATCH_DIR}/${BATCH_NAME}.err" ]] && [[ -s "${BATCH_DIR}/${BATCH_NAME}.err" ]]; then
                    FAILED=$((FAILED + 1))
                    echo "  [FAILED]  ${BATCH_NAME} — check ${BATCH_DIR}/${BATCH_NAME}.err"
                else
                    DONE=$((DONE + 1))
                    echo "  [DONE?]   ${BATCH_NAME} (no .done marker, job completed)"
                fi
            fi
        else
            echo "  [PENDING] ${BATCH_NAME}"
        fi
    done

    echo ""
    echo "Total: ${TOTAL} | Done: ${DONE} | Running: ${RUNNING} | Failed: ${FAILED}"

    if [[ $DONE -eq $TOTAL ]]; then
        echo ""
        echo "All batches complete! Run: $0 --concat"
    fi
fi

# --- Concatenate results ---
if [[ "$ACTION" == "concat" ]]; then
    echo "=== Concatenating PRS Results ==="

    SCORE_FILES=()
    # pgsc_calc v2.2.0 outputs to: batch_XXX/SAMPLENAME/score/aggregated_scores.txt.gz
    while IFS= read -r -d '' f; do
        SCORE_FILES+=("$f")
    done < <(find "${RESULTS_DIR}" -name "aggregated_scores.txt.gz" -print0 2>/dev/null)

    if [[ ${#SCORE_FILES[@]} -eq 0 ]]; then
        echo "No score files found in ${RESULTS_DIR}/"
        echo "Expected pgsc_calc output files (aggregated_scores.txt.gz)"
        exit 1
    fi

    echo "Found ${#SCORE_FILES[@]} score files"

    OUTPUT="${RESULTS_DIR}/all_scores.txt.gz"

    # Concatenate: keep header from first file, skip headers from rest
    {
        zcat "${SCORE_FILES[0]}" | head -1  # Header
        for f in "${SCORE_FILES[@]}"; do
            zcat "$f" | tail -n +2  # Data without header
        done
    } | gzip > "$OUTPUT"

    TOTAL_LINES=$(($(zcat "$OUTPUT" | wc -l) - 1))
    echo "Concatenated ${TOTAL_LINES} score lines to ${OUTPUT}"
    echo ""
    echo "Next: Copy to local machine and run process_prs_results.py"
    echo "  scp YOUR_HPC_HOST:\${OUTPUT} /path/to/local/GenomeAnalysis/results/"
fi
