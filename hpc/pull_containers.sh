#!/usr/bin/env bash
set -euo pipefail

# Pull BioContainers for triage tools (run on login node — needs internet)
# Usage: nohup bash pull_containers.sh > pull_containers.log 2>&1 &

CONTAINER_DIR="/rds/general/user/${HPC_USER:-YOUR_USERNAME}/home/GenomeAnalysis/containers"
mkdir -p "${CONTAINER_DIR}"

echo "=== Pulling triage containers ==="
echo "Output: ${CONTAINER_DIR}"
echo "Start: $(date)"
echo

# vcfanno
if [[ ! -f "${CONTAINER_DIR}/vcfanno_0.3.7.sif" ]]; then
    echo "[1/3] Pulling vcfanno 0.3.7..."
    apptainer pull "${CONTAINER_DIR}/vcfanno_0.3.7.sif" \
        docker://quay.io/biocontainers/vcfanno:0.3.7--he881be0_0
    echo "  -> Done"
else
    echo "[1/3] vcfanno already pulled"
fi
echo

# slivar
if [[ ! -f "${CONTAINER_DIR}/slivar_0.3.2.sif" ]]; then
    echo "[2/3] Pulling slivar 0.3.2..."
    apptainer pull "${CONTAINER_DIR}/slivar_0.3.2.sif" \
        docker://quay.io/biocontainers/slivar:0.3.2--h5f107b1_0
    echo "  -> Done"
else
    echo "[2/3] slivar already pulled"
fi
echo

# VEP (if not already pulled by run_vep.sh)
if [[ ! -f "${CONTAINER_DIR}/ensembl-vep_112.0.sif" ]]; then
    echo "[3/3] Pulling ensembl-vep 112.0..."
    apptainer pull "${CONTAINER_DIR}/ensembl-vep_112.0.sif" \
        docker://ensemblorg/ensembl-vep:release_112.0
    echo "  -> Done"
else
    echo "[3/3] ensembl-vep already pulled"
fi
echo

echo "=== All containers pulled ==="
echo "End: $(date)"
ls -lh "${CONTAINER_DIR}"/*.sif
