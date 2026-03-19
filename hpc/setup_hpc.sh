#!/usr/bin/env bash
set -euo pipefail

# One-time HPC environment setup for nf-core/raredisease pipeline
# Imperial College HPC (CX3) — PBS Pro + Apptainer
#
# Usage: ssh imperial 'bash -s' < setup_hpc.sh

echo "=== Imperial HPC Setup for GenomeAnalysis ==="
echo "User: $(whoami)"
echo "Host: $(hostname)"
echo "Date: $(date)"
echo ""

# ─── 0. Initialize Module System ─────────────────────────────
echo "--- Initializing module system ---"

source /usr/share/lmod/lmod/init/bash 2>/dev/null
export MODULEPATH=/rds/general/apps/modules/modulegroups:/rds/easybuild/noarch/apps/modules/all:/etc/modulefiles:/usr/share/modulefiles
module load tools/prod 2>/dev/null || true
module load tools/bioinf 2>/dev/null || true

if type module &>/dev/null; then
    echo "  Lmod initialized successfully"
else
    echo "  WARNING: Module system not available"
fi
echo ""

# ─── 1. Directory Structure ─────────────────────────────────
echo "--- Creating directory structure ---"

DIRS=(
    "/rds/general/user/${HPC_USER:-YOUR_USERNAME}/ephemeral/singularity_cache"
    "/rds/general/user/${HPC_USER:-YOUR_USERNAME}/ephemeral/nf-work"
    "/rds/general/user/${HPC_USER:-YOUR_USERNAME}/ephemeral/fastq"
    "/rds/general/user/${HPC_USER:-YOUR_USERNAME}/ephemeral/references"
    "/rds/general/user/${HPC_USER:-YOUR_USERNAME}/home/GenomeAnalysis/results/pipeline_info"
)

for d in "${DIRS[@]}"; do
    if [[ ! -d "${d}" ]]; then
        mkdir -p "${d}"
        echo "  Created: ${d}"
    else
        echo "  Exists:  ${d}"
    fi
done
echo ""

# ─── 2. Check Available Modules ──────────────────────────────
echo "--- Checking available modules ---"

MODULES_TO_CHECK=(
    "Java/17"
    "Nextflow/25.10.2"
    "SAMtools/1.21"
    "BCFtools/1.21"
    "BEDTools/2.31"
)

for mod in "${MODULES_TO_CHECK[@]}"; do
    if module -t avail "${mod}" 2>&1 | grep -q "${mod}"; then
        echo "  Available: ${mod}"
    else
        echo "  Not found: ${mod} (may have different name)"
    fi
done
echo ""

# ─── 3. Install/Load Nextflow ─────────────────────────────────
echo "--- Setting up Nextflow ---"

module load Java/17 2>/dev/null || module load Java/17.0.6 2>/dev/null || true
module load Nextflow/25.10.2 2>/dev/null || module load Nextflow/24.04.2 2>/dev/null || true

if command -v nextflow &>/dev/null; then
    echo "  Nextflow available via module: $(nextflow -version 2>&1 | grep 'version' | head -1)"
else
    echo "  Nextflow module not loaded, installing to ~/bin..."
    NF_DIR="${HOME}/bin"
    mkdir -p "${NF_DIR}"
    if [[ -f "${NF_DIR}/nextflow" ]]; then
        echo "  Binary exists at ${NF_DIR}/nextflow"
    else
        cd "${NF_DIR}"
        curl -s https://get.nextflow.io | bash
        chmod +x nextflow
        echo "  Installed: ${NF_DIR}/nextflow"
    fi
    export PATH="${NF_DIR}:${PATH}"
fi
echo ""

# ─── 4. Configure Environment in .bashrc ─────────────────────
echo "--- Configuring environment ---"

BASHRC="${HOME}/.bashrc"

# Module initialization for non-interactive shells
if ! grep -q "lmod/init/bash" "${BASHRC}" 2>/dev/null; then
    cat >> "${BASHRC}" << 'LMOD_INIT'

# === GenomeAnalysis pipeline environment ===
# Initialize Lmod for non-interactive shells
source /usr/share/lmod/lmod/init/bash 2>/dev/null
export MODULEPATH=/rds/general/apps/modules/modulegroups:/rds/easybuild/noarch/apps/modules/all:/etc/modulefiles:/usr/share/modulefiles
module load tools/prod 2>/dev/null
module load tools/bioinf 2>/dev/null
module load Java/17 2>/dev/null
module load Nextflow/25.10.2 2>/dev/null || module load Nextflow/24.04.2 2>/dev/null
LMOD_INIT
    echo "  Added Lmod initialization to ~/.bashrc"
else
    echo "  Lmod init already in ~/.bashrc"
fi

# Singularity/Apptainer cache
SINGULARITY_CACHEDIR="/rds/general/user/${HPC_USER:-YOUR_USERNAME}/ephemeral/singularity_cache"
if ! grep -q "SINGULARITY_CACHEDIR" "${BASHRC}" 2>/dev/null; then
    cat >> "${BASHRC}" << EOF

# Singularity/Apptainer cache on ephemeral (large containers)
export SINGULARITY_CACHEDIR="${SINGULARITY_CACHEDIR}"
export NXF_SINGULARITY_CACHEDIR="${SINGULARITY_CACHEDIR}"
export APPTAINER_CACHEDIR="${SINGULARITY_CACHEDIR}"
EOF
    echo "  Added SINGULARITY_CACHEDIR to ~/.bashrc"
else
    echo "  SINGULARITY_CACHEDIR already in ~/.bashrc"
fi

# Nextflow work directory
if ! grep -q "NXF_WORK" "${BASHRC}" 2>/dev/null; then
    cat >> "${BASHRC}" << 'EOF'

# Nextflow work directory on ephemeral
export NXF_WORK="/rds/general/user/${HPC_USER:-YOUR_USERNAME}/ephemeral/nf-work"
EOF
    echo "  Added NXF_WORK to ~/.bashrc"
else
    echo "  NXF_WORK already in ~/.bashrc"
fi
echo ""

# ─── 5. Test PBS Submission ───────────────────────────────────
echo "--- Testing PBS submission ---"

if command -v qsub &>/dev/null || [[ -x /opt/pbs/bin/qsub ]]; then
    QSUB="${QSUB:-$(command -v qsub 2>/dev/null || echo /opt/pbs/bin/qsub)}"
    JOBID=$("${QSUB}" -N nf-test -l walltime=00:05:00,select=1:ncpus=1:mem=1gb \
        -o /tmp/nf-test.out -e /tmp/nf-test.err \
        -- /bin/bash -c "echo 'Hello from PBS'; hostname; date; singularity --version 2>/dev/null; java -version 2>&1 | head -1" 2>&1)
    echo "  PBS test job submitted: ${JOBID}"
    echo "  Check with: qstat -u ${HPC_USER:-YOUR_USERNAME}"
else
    echo "  qsub not available"
fi
echo ""

# ─── 6. Test Singularity/Apptainer ───────────────────────────
echo "--- Testing Singularity/Apptainer ---"

if command -v singularity &>/dev/null; then
    echo "  Singularity/Apptainer version: $(singularity --version)"
    echo "  Cache dir: ${SINGULARITY_CACHEDIR}"
else
    echo "  Singularity not found in PATH"
fi
echo ""

# ─── 7. Verify Nextflow ──────────────────────────────────────
echo "--- Verifying Nextflow ---"

if command -v nextflow &>/dev/null; then
    nextflow -version 2>&1 | head -3
else
    echo "  Nextflow not in PATH"
fi
echo ""

# ─── 8. Storage Check ────────────────────────────────────────
echo "=== Storage Summary ==="
echo ""
echo "Home:"
quota -s 2>/dev/null || df -h /rds/general/user/${HPC_USER:-YOUR_USERNAME}/home/ 2>/dev/null || echo "Cannot check home quota"
echo ""
echo "Ephemeral:"
df -h /rds/general/user/${HPC_USER:-YOUR_USERNAME}/ephemeral/ 2>/dev/null || echo "Cannot check ephemeral storage"
echo ""

# ─── Summary ─────────────────────────────────────────────────
echo "=== Setup Complete ==="
echo ""
echo "Scheduler: PBS Pro (qsub/qstat)"
echo "Container: Apptainer (singularity CLI compatible)"
echo "Modules: tools/prod + tools/bioinf loaded"
echo ""
echo "Directory structure:"
for d in "${DIRS[@]}"; do
    echo "  ${d}"
done
echo ""
echo "Next steps:"
echo "  1. source ~/.bashrc"
echo "  2. Run stage_references.sh to download reference data"
echo "  3. Run download_fastq.sh to get FASTQ files"
echo "  4. Run run_pipeline.sh to launch the pipeline"
echo ""
echo "Estimated storage needs:"
echo "  References:  ~120-150 GB"
echo "  FASTQs:      ~100-150 GB"
echo "  Work dir:    ~200-500 GB (temporary, cleaned on success)"
echo "  Results:     ~50-100 GB"
echo "  Containers:  ~20-30 GB"
echo "  TOTAL:       ~500 GB - 1 TB ephemeral"
