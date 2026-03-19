#!/usr/bin/env bash
set -euo pipefail

# Download triage reference files for GRCh38 variant annotation
# Run on HPC login node (compute nodes have no internet)
#
# Usage: nohup bash download_triage_refs.sh > download_triage_refs.log 2>&1 &

source /usr/share/lmod/lmod/init/bash
module purge 2>/dev/null || true
module load tools/prod SAMtools/1.21-GCC-13.3.0 HTSlib/1.21-GCC-13.3.0

# Verify bgzip/tabix are available (from HTSlib, not SAMtools)
command -v bgzip &>/dev/null || { echo "FATAL: bgzip not found"; exit 1; }
echo "Using bgzip: $(which bgzip)"
echo "Using tabix: $(which tabix)"

REFDIR="/rds/general/user/${HPC_USER:-YOUR_USERNAME}/home/GenomeAnalysis/references/triage"
mkdir -p "${REFDIR}"

echo "=== Downloading triage reference files (GRCh38) ==="
echo "Output: ${REFDIR}"
echo "Start: $(date)"
echo ""

# ============================================================================
# 1. Segmental duplications (UCSC genomicSuperDups)
# ============================================================================
echo "[1/5] Segmental duplications..."
if [[ ! -f "${REFDIR}/genomicSuperDups.bed.gz" ]]; then
    wget -q -O - "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz" \
        | gunzip -c \
        | awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$5}' \
        | sort -k1,1 -k2,2n \
        | bgzip -c > "${REFDIR}/genomicSuperDups.bed.gz"
    tabix -p bed "${REFDIR}/genomicSuperDups.bed.gz"
    echo "  -> Done: $(wc -l < <(zcat "${REFDIR}/genomicSuperDups.bed.gz")) regions"
else
    echo "  -> Already exists, skipping"
fi
echo ""

# ============================================================================
# 2. Low-complexity regions (Heng Li's mdust LCR track for GRCh38)
# ============================================================================
echo "[2/5] Low-complexity regions..."
if [[ ! -f "${REFDIR}/LCR-hs38.bed.gz" ]]; then
    wget -q -O - "https://github.com/lh3/varcmp/raw/master/scripts/LCR-hs38.bed.gz" \
        > "${REFDIR}/LCR-hs38.bed.gz"
    tabix -p bed "${REFDIR}/LCR-hs38.bed.gz"
    echo "  -> Done"
else
    echo "  -> Already exists, skipping"
fi
echo ""

# ============================================================================
# 3. Simple repeats (UCSC simpleRepeat track)
# ============================================================================
echo "[3/5] Simple repeats..."
if [[ ! -f "${REFDIR}/simpleRepeat.bed.gz" ]]; then
    wget -q -O - "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz" \
        | gunzip -c \
        | awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$17}' \
        | sort -k1,1 -k2,2n \
        | bgzip -c > "${REFDIR}/simpleRepeat.bed.gz"
    tabix -p bed "${REFDIR}/simpleRepeat.bed.gz"
    echo "  -> Done"
else
    echo "  -> Already exists, skipping"
fi
echo ""

# ============================================================================
# 4. Tandem repeat finder (UCSC TRF — same as simpleRepeat but different format)
# ============================================================================
echo "[4/5] Tandem repeat finder regions..."
if [[ ! -f "${REFDIR}/trf.bed.gz" ]]; then
    # TRF track is embedded in simpleRepeat on UCSC for hg38
    # Extract just chr, start, end, period for flagging
    wget -q -O - "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz" \
        | gunzip -c \
        | awk 'BEGIN{OFS="\t"} {print $2,$3,$4,"TRF_"$6}' \
        | sort -k1,1 -k2,2n \
        | bgzip -c > "${REFDIR}/trf.bed.gz"
    tabix -p bed "${REFDIR}/trf.bed.gz"
    echo "  -> Done"
else
    echo "  -> Already exists, skipping"
fi
echo ""

# ============================================================================
# 5. gnomAD v4.1 genome sites (AF-only VCF — much smaller than full sites)
# ============================================================================
echo "[5/5] gnomAD v4.1 AF-only sites..."
echo "  NOTE: This downloads per-chromosome VCFs and merges them."
echo "  Full gnomAD sites is ~700GB. We download the AF-only version (~60GB total)."
echo ""

GNOMAD_DIR="${REFDIR}/gnomad"
mkdir -p "${GNOMAD_DIR}"

# gnomAD v4.1 provides joint frequency data
# Using the sites VCF with only AF fields to keep size manageable
GNOMAD_BASE="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes"

FAILED=0
for chr in $(seq 1 22) X Y; do
    OUTFILE="${GNOMAD_DIR}/gnomad.genomes.v4.1.sites.chr${chr}.vcf.bgz"
    if [[ -f "${OUTFILE}" && -f "${OUTFILE}.tbi" ]]; then
        echo "  chr${chr}: already downloaded"
        continue
    fi
    echo "  chr${chr}: downloading..."
    URL="${GNOMAD_BASE}/gnomad.genomes.v4.1.sites.chr${chr}.vcf.bgz"
    if wget -q -c -O "${OUTFILE}" "${URL}"; then
        # Download index too
        wget -q -c -O "${OUTFILE}.tbi" "${URL}.tbi" || true
        echo "  chr${chr}: done ($(ls -lh "${OUTFILE}" | awk '{print $5}'))"
    else
        echo "  chr${chr}: FAILED"
        FAILED=$((FAILED + 1))
        rm -f "${OUTFILE}"
    fi
done

if [[ ${FAILED} -gt 0 ]]; then
    echo ""
    echo "WARNING: ${FAILED} chromosome(s) failed to download."
    echo "Re-run this script to retry (wget -c resumes partial downloads)."
fi

echo ""
echo "=== Download complete ==="
echo "End: $(date)"
echo ""
echo "Files in ${REFDIR}:"
ls -lh "${REFDIR}"/*.bed.gz 2>/dev/null || true
echo ""
echo "gnomAD files:"
ls -lh "${GNOMAD_DIR}"/*.vcf.bgz 2>/dev/null | head -5
echo "... ($(ls "${GNOMAD_DIR}"/*.vcf.bgz 2>/dev/null | wc -l) files total)"
echo ""
echo "Next steps:"
echo "  1. Update vcfanno.toml paths to point to ${REFDIR}/"
echo "  2. For gnomAD, either merge chromosomes or update vcfanno to use per-chr"
echo "  3. Install vcfanno and slivar (check if available as modules or use containers)"
