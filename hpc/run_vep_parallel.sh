#!/bin/bash
#
# Submit per-chromosome VEP annotation jobs in parallel
#
# Usage: bash run_vep_parallel.sh [INPUT_VCF]
#
# Submits 24 PBS jobs (chr1-22, chrX, chrY) each running VEP on one chromosome.
# Each job requests 4h walltime (vs 24h+ for whole genome).
# After all jobs complete, run merge_vep.sh to combine results.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$(readlink -f /rds/general/user/nskene/home/GenomeAnalysis)"
RESULTS="${BASE}/results"

INPUT_VCF="$(readlink -f "${1:-${RESULTS}/call_snv/genome/GFXC087577_case_snv.vcf.gz}")"

if [[ ! -f "${INPUT_VCF}" ]]; then
    echo "ERROR: Cannot find DeepVariant VCF at: ${INPUT_VCF}"
    echo "Usage: bash run_vep_parallel.sh [/path/to/input.vcf.gz]"
    exit 1
fi

# Chromosomes to annotate
CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)

# Output and log directories
mkdir -p "${RESULTS}/vep/per_chr" "${RESULTS}/vep/logs"

echo "============================================================"
echo "VEP Parallel Submission"
echo "============================================================"
echo "Input:  ${INPUT_VCF}"
echo "Jobs:   ${#CHROMS[@]} chromosomes"
echo "Date:   $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
echo

# Check which chromosomes already have completed VEP output
SKIP=0
SUBMIT=0
JOB_IDS=()

for CHR in "${CHROMS[@]}"; do
    OUTPUT="${RESULTS}/vep/per_chr/GFXC087577.vep.${CHR}.vcf.gz"

    # Skip if already annotated (has index = completed successfully)
    if [[ -f "${OUTPUT}" && -f "${OUTPUT}.tbi" ]]; then
        echo "  SKIP  ${CHR}  (already annotated: ${OUTPUT})"
        SKIP=$((SKIP + 1))
        continue
    fi

    # Submit PBS job
    JOB_ID=$(qsub \
        -N "vep_${CHR}" \
        -v "CHR=${CHR},INPUT_VCF=${INPUT_VCF}" \
        -o "${RESULTS}/vep/logs/vep_${CHR}.out" \
        -e "${RESULTS}/vep/logs/vep_${CHR}.err" \
        "${SCRIPT_DIR}/run_vep.sh")

    echo "  SUBMIT  ${CHR}  -> ${JOB_ID}"
    JOB_IDS+=("${JOB_ID}")
    SUBMIT=$((SUBMIT + 1))
done

echo
echo "============================================================"
echo "Summary"
echo "============================================================"
echo "Submitted: ${SUBMIT} jobs"
echo "Skipped:   ${SKIP} (already complete)"
echo "Total:     ${#CHROMS[@]} chromosomes"
echo

if [[ ${SUBMIT} -gt 0 ]]; then
    echo "Monitor progress:"
    echo "  qstat -u nskene | grep vep_"
    echo
    echo "Check for failures:"
    echo "  grep -l ERROR ${RESULTS}/vep/logs/vep_*.err"
    echo
    echo "After all jobs complete, merge results:"
    echo "  bash ${SCRIPT_DIR}/merge_vep.sh"
fi
