#!/bin/bash
#
# run_pgscalc.sh — Launch pgsc_calc on Imperial HPC (CX3)
#
# Single-batch run of pgsc_calc Nextflow pipeline.
# For multi-batch runs across all PGS Catalog scores, use prs_batch.sh instead.
#
# Usage:
#   bash run_pgscalc.sh PGS000001,PGS000002,PGS000003
#   bash run_pgscalc.sh --file pgs_ids_batch1.txt
#   bash run_pgscalc.sh --file pgs_ids_batch1.txt --resume
#   bash run_pgscalc.sh --file pgs_ids_batch1.txt --batch-name batch_000
#

set -euo pipefail

# --- Configuration ---
HPC_BASE="/rds/general/user/${HPC_USER:-YOUR_USERNAME}/home/GenomeAnalysis"
WORK_BASE="/rds/general/user/${HPC_USER:-YOUR_USERNAME}/ephemeral/nf-work-prs"
CACHE_DIR="/rds/general/user/${HPC_USER:-YOUR_USERNAME}/ephemeral/singularity_cache"
RESULTS_DIR="${HPC_BASE}/results/prs"
SCRIPTS_DIR="${HPC_BASE}/scripts"
SAMPLESHEET="${SCRIPTS_DIR}/prs_samplesheet.csv"

# pgsc_calc version
PGSC_CALC_VERSION="v2.2.0"

# Target build — matches current VCF (GRCh37)
TARGET_BUILD="GRCh37"

# --- Parse arguments ---
PGS_IDS=""
PGS_FILE=""
RESUME=""
BATCH_NAME=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --file|-f)
            PGS_FILE="$2"
            shift 2
            ;;
        --resume|-r)
            RESUME="-resume"
            shift
            ;;
        --build)
            TARGET_BUILD="$2"
            shift 2
            ;;
        --batch-name)
            BATCH_NAME="$2"
            shift 2
            ;;
        *)
            PGS_IDS="$1"
            shift
            ;;
    esac
done

# Read PGS IDs from file if provided
if [[ -n "$PGS_FILE" ]]; then
    if [[ ! -f "$PGS_FILE" ]]; then
        echo "ERROR: PGS ID file not found: $PGS_FILE" >&2
        exit 1
    fi
    # Read file, join with commas, strip whitespace
    PGS_IDS=$(tr '\n' ',' < "$PGS_FILE" | sed 's/,$//' | sed 's/ //g')
    # Derive batch name from file if not specified
    if [[ -z "$BATCH_NAME" ]]; then
        BATCH_NAME=$(basename "$PGS_FILE" .txt)
    fi
fi

if [[ -z "$PGS_IDS" ]]; then
    echo "Usage: $0 PGS000001,PGS000002 [--resume]"
    echo "       $0 --file pgs_ids.txt [--resume] [--batch-name NAME]"
    exit 1
fi

# Per-batch work and output directories to avoid lock contention
WORK_DIR="${WORK_BASE}/${BATCH_NAME:-default}"
OUTDIR="${RESULTS_DIR}/${BATCH_NAME:-default}"

# Count IDs
N_IDS=$(echo "$PGS_IDS" | tr ',' '\n' | wc -l | tr -d ' ')
echo "=== pgsc_calc Launch ==="
echo "PGS IDs: ${N_IDS}"
echo "Batch: ${BATCH_NAME:-default}"
echo "Target build: ${TARGET_BUILD}"
echo "Results: ${OUTDIR}"
echo "Work dir: ${WORK_DIR}"
echo ""

# --- Environment setup ---
source /usr/share/lmod/lmod/init/bash
module purge
module load tools/prod
module load Java/17
module load Nextflow/25.10.2

# Apptainer/Singularity cache (shared across batches is fine)
export NXF_SINGULARITY_CACHEDIR="${CACHE_DIR}"
export SINGULARITY_CACHEDIR="${CACHE_DIR}"

# Per-batch NXF_HOME to avoid .nextflow/cache lock contention between parallel batches
export NXF_HOME="${WORK_DIR}/.nextflow-home"

# Create directories
mkdir -p "${OUTDIR}" "${WORK_DIR}" "${CACHE_DIR}" "${NXF_HOME}"

# --- Launch pgsc_calc ---
# No custom nextflow.config — pgsc_calc has its own config.
# Our raredisease nextflow.config has params that conflict with pgsc_calc.
echo "Launching pgsc_calc ${PGSC_CALC_VERSION}..."
echo "Start time: $(date)"
echo ""

nextflow run pgscatalog/pgsc_calc \
    -r "${PGSC_CALC_VERSION}" \
    -profile singularity \
    --input "${SAMPLESHEET}" \
    --pgs_id "${PGS_IDS}" \
    --target_build "${TARGET_BUILD}" \
    --outdir "${OUTDIR}" \
    -work-dir "${WORK_DIR}" \
    ${RESUME}

echo ""
echo "pgsc_calc finished: $(date)"
echo "Results in: ${OUTDIR}"
