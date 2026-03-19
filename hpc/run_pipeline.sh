#!/usr/bin/env bash
set -euo pipefail

# Launch nf-core/raredisease v2.6.0 on Imperial HPC (CX3, PBS Pro)
# Sample: ${SAMPLE_ID:-YOUR_SAMPLE} (Dante Labs WGS, singleton)
#
# Module initialization:
#   source /usr/share/lmod/lmod/init/bash
#   module purge
#   module load tools/prod Java/17 Nextflow/25.10.2

# ─── Configuration ───────────────────────────────────────────
PIPELINE="nf-core/raredisease"
REVISION="2.6.0"
PROFILE="imperial"

# Ensure PBS commands (qsub, qstat) are in PATH (login-bi doesn't have them)
export PATH="/opt/pbs/bin:${PATH}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${SCRIPT_DIR}/nextflow.config"
SAMPLESHEET="${SCRIPT_DIR}/samplesheet.csv"

REFDIR="/rds/general/user/${HPC_USER:-YOUR_USERNAME}/home/GenomeAnalysis/references"
OUTDIR="/rds/general/user/${HPC_USER:-YOUR_USERNAME}/home/GenomeAnalysis/results"
WORKDIR="/rds/general/user/${HPC_USER:-YOUR_USERNAME}/ephemeral/nf-work"

# Isolate Nextflow metadata from other pipelines (e.g. PRS batches)
LAUNCH_DIR="/rds/general/user/${HPC_USER:-YOUR_USERNAME}/home/GenomeAnalysis/raredisease_run"
mkdir -p "${LAUNCH_DIR}"
export NXF_HOME="${LAUNCH_DIR}/.nextflow-home"

# ─── Pre-flight checks ──────────────────────────────────────
echo "=== nf-core/raredisease v${REVISION} Launch ==="
echo ""

# Check Nextflow is available
if ! command -v nextflow &>/dev/null; then
    echo "ERROR: nextflow not found. Run setup_hpc.sh first."
    exit 1
fi

# Check config and samplesheet exist
for f in "${CONFIG}" "${SAMPLESHEET}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: Missing file: ${f}"
        exit 1
    fi
done

# Check reference files
echo "Checking reference data..."
MISSING=0
for f in \
    "${REFDIR}/genome/Homo_sapiens_assembly38.fasta" \
    "${REFDIR}/vep/homo_sapiens" \
    "${REFDIR}/cadd/whole_genome_SNVs.tsv.gz" \
; do
    if [[ ! -e "${f}" ]]; then
        echo "  MISSING: ${f}"
        MISSING=$((MISSING + 1))
    fi
done

if [[ ${MISSING} -gt 0 ]]; then
    echo "WARNING: ${MISSING} reference file(s) missing. Run stage_references.sh first."
    if [[ -t 0 ]]; then
        echo "Continue anyway? (y/N)"
        read -r response
        if [[ "${response}" != "y" && "${response}" != "Y" ]]; then
            exit 1
        fi
    else
        echo "Non-interactive mode: aborting due to missing references."
        exit 1
    fi
fi

# Check FASTQs exist
if [[ ! -f "/rds/general/user/${HPC_USER:-YOUR_USERNAME}/home/GenomeAnalysis/fastq/${SAMPLE_ID:-YOUR_SAMPLE}_R1.fastq.gz" ]]; then
    echo "ERROR: FASTQ files not found. Run download_fastq.sh first."
    exit 1
fi

echo "Pre-flight checks passed."
echo ""

# ─── Create output directories ──────────────────────────────
mkdir -p "${OUTDIR}/pipeline_info" "${WORKDIR}"

# ─── Launch pipeline ────────────────────────────────────────
echo "Launching pipeline..."
echo "  Pipeline:    ${PIPELINE} -r ${REVISION}"
echo "  Profile:     ${PROFILE}"
echo "  Samplesheet: ${SAMPLESHEET}"
echo "  Output:      ${OUTDIR}"
echo "  Work dir:    ${WORKDIR}"
echo "  Launch dir:  ${LAUNCH_DIR}"
echo ""

# Run from isolated launch dir so .nextflow/ cache doesn't collide with PRS batches
cd "${LAUNCH_DIR}"

nextflow run "${PIPELINE}" \
    -r "${REVISION}" \
    -profile "${PROFILE}" \
    -c "${CONFIG}" \
    --input "${SAMPLESHEET}" \
    --genome GRCh38 \
    --outdir "${OUTDIR}" \
    -work-dir "${WORKDIR}" \
    \
    --analysis_type wgs \
    --variant_caller deepvariant \
    --fasta "${REFDIR}/genome/Homo_sapiens_assembly38.fasta" \
    --intervals_wgs "/rds/general/user/${HPC_USER:-YOUR_USERNAME}/home/.nextflow/assets/nf-core/raredisease/assets/chr1-chr22chrXchrYchrM_grch38.interval_list" \
    --intervals_y "/rds/general/user/${HPC_USER:-YOUR_USERNAME}/home/.nextflow/assets/nf-core/raredisease/assets/chrY_grch38.interval_list" \
    --vep_cache "${REFDIR}/vep" \
    --vep_cache_version 114 \
    --cadd_resources "${REFDIR}/cadd" \
    --skip_subworkflows me_calling,me_annotation,mt_annotation,repeat_annotation,repeat_calling,sv_annotation,sv_calling,snv_annotation,generate_clinical_set \
    --skip_tools gens,germlinecnvcaller,vcf2cytosure,eklipse \
    \
    -resume \
    \
    2>&1 | tee "${OUTDIR}/pipeline_info/nextflow_stdout.log"

EXIT_CODE=$?

echo ""
if [[ ${EXIT_CODE} -eq 0 ]]; then
    echo "=== Pipeline completed successfully ==="
    echo "Results: ${OUTDIR}"
    echo "Reports: ${OUTDIR}/pipeline_info/"
else
    echo "=== Pipeline failed with exit code ${EXIT_CODE} ==="
    echo "Check logs:"
    echo "  Nextflow log: .nextflow.log"
    echo "  Stdout log:   ${OUTDIR}/pipeline_info/nextflow_stdout.log"
    echo "  Trace:        ${OUTDIR}/pipeline_info/trace.txt"
    echo ""
    echo "To resume from last checkpoint: re-run this script (uses -resume)"
fi

exit ${EXIT_CODE}
