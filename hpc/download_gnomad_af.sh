#!/usr/bin/env bash
set -euo pipefail

# Download gnomAD v4.1 genome sites and extract AF-only fields
# Downloads one chromosome at a time, strips all INFO except AF fields, deletes full file
#
# Result: ~30-60 GB AF-only VCFs (vs ~700 GB full sites)
# Peak temp: ~45 GB (one full chr + partial output)
#
# Usage: nohup bash download_gnomad_af.sh > download_gnomad_af.log 2>&1 &
# Run on HPC login node (compute nodes have no internet)

source /usr/share/lmod/lmod/init/bash
module purge 2>/dev/null || true
module load tools/prod BCFtools/1.21-GCC-13.3.0 HTSlib/1.21-GCC-13.3.0

command -v bcftools &>/dev/null || { echo "FATAL: bcftools not found"; exit 1; }
command -v bgzip &>/dev/null || { echo "FATAL: bgzip not found"; exit 1; }
command -v tabix &>/dev/null || { echo "FATAL: tabix not found"; exit 1; }

OUTDIR="/rds/general/user/${HPC_USER:-YOUR_USERNAME}/home/GenomeAnalysis/references/triage/gnomad"
TMPDIR_DL="${OUTDIR}/tmp_download"
mkdir -p "${OUTDIR}" "${TMPDIR_DL}"

GNOMAD_BASE="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes"

# Fields to keep from gnomAD v4.1:
#   AF          - overall allele frequency
#   AF_grpmax   - max AF across genetic ancestry groups (replaces AF_popmax in v3)
#   AC          - allele count
#   AN          - allele number
#   nhomalt     - number of homozygous alternate individuals
#   fafmax_faf95_max - max filtering AF at 95% CI across groups (recommended for ACMG BA1/BS1)
KEEP_FIELDS="INFO/AF,INFO/AF_grpmax,INFO/AC,INFO/AN,INFO/nhomalt,INFO/fafmax_faf95_max"

echo "=== gnomAD v4.1 AF-Only Extraction ==="
echo "Output: ${OUTDIR}"
echo "Keeping: ${KEEP_FIELDS}"
echo "Start: $(date)"
echo ""

# First: verify field names exist by downloading chr22 header (smallest chr)
echo "[Pre-flight] Checking gnomAD v4.1 field names on chr22..."
HEADER_URL="${GNOMAD_BASE}/gnomad.genomes.v4.1.sites.chr22.vcf.bgz"
HEADER_FILE="${TMPDIR_DL}/chr22_header_check.vcf.bgz"

if [[ ! -f "${OUTDIR}/gnomad.v4.1.af_only.chr22.vcf.bgz" ]]; then
    # Download just the first few bytes to get the header
    # bcftools can read headers from remote URLs if htslib supports it
    # But safer to download the tbi first, then use bcftools view -h with region to get header only
    echo "  Downloading chr22 index for header check..."
    wget -q -O "${TMPDIR_DL}/chr22.vcf.bgz.tbi" "${HEADER_URL}.tbi" 2>/dev/null || true

    echo "  Downloading chr22 (smallest, ~10GB)..."
    wget -q -c -O "${HEADER_FILE}" "${HEADER_URL}"

    echo "  Checking INFO fields in header..."
    FIELDS_FOUND=0
    for field in AF AF_grpmax AC AN nhomalt fafmax_faf95_max; do
        if bcftools view -h "${HEADER_FILE}" | grep -q "ID=${field},"; then
            echo "    ${field}: FOUND"
            FIELDS_FOUND=$((FIELDS_FOUND + 1))
        else
            echo "    ${field}: NOT FOUND"
        fi
    done

    if [[ ${FIELDS_FOUND} -lt 4 ]]; then
        echo ""
        echo "WARNING: Some expected fields not found in gnomAD v4.1 header."
        echo "Full INFO fields available:"
        bcftools view -h "${HEADER_FILE}" | grep "^##INFO=" | sed 's/.*ID=//; s/,.*//' | tr '\n' ' '
        echo ""
        echo "Proceeding with available fields..."
    fi
    echo ""

    # Process chr22 as first chromosome
    echo "[chr22] Extracting AF-only fields..."
    bcftools annotate \
        -x "^${KEEP_FIELDS}" \
        -Oz -o "${OUTDIR}/gnomad.v4.1.af_only.chr22.vcf.bgz" \
        "${HEADER_FILE}"
    tabix -p vcf "${OUTDIR}/gnomad.v4.1.af_only.chr22.vcf.bgz"

    FULL_SIZE=$(ls -lh "${HEADER_FILE}" | awk '{print $5}')
    AF_SIZE=$(ls -lh "${OUTDIR}/gnomad.v4.1.af_only.chr22.vcf.bgz" | awk '{print $5}')
    echo "  Full: ${FULL_SIZE} -> AF-only: ${AF_SIZE}"

    # Clean up full file
    rm -f "${HEADER_FILE}" "${TMPDIR_DL}/chr22.vcf.bgz.tbi"
    echo "  Deleted full sites VCF"
    echo ""
fi

# Process remaining chromosomes
FAILED=0
for chr in $(seq 1 21) X Y; do
    AF_FILE="${OUTDIR}/gnomad.v4.1.af_only.chr${chr}.vcf.bgz"

    if [[ -f "${AF_FILE}" && -f "${AF_FILE}.tbi" ]]; then
        echo "[chr${chr}] Already done, skipping"
        continue
    fi

    echo "[chr${chr}] Downloading full sites VCF..."
    FULL_FILE="${TMPDIR_DL}/gnomad.genomes.v4.1.sites.chr${chr}.vcf.bgz"

    if ! wget -q -c -O "${FULL_FILE}" "${GNOMAD_BASE}/gnomad.genomes.v4.1.sites.chr${chr}.vcf.bgz"; then
        echo "  FAILED to download chr${chr}"
        FAILED=$((FAILED + 1))
        rm -f "${FULL_FILE}"
        continue
    fi

    FULL_SIZE=$(ls -lh "${FULL_FILE}" | awk '{print $5}')
    echo "  Downloaded: ${FULL_SIZE}"

    echo "  Extracting AF-only fields..."
    if bcftools annotate \
        -x "^${KEEP_FIELDS}" \
        -Oz -o "${AF_FILE}" \
        "${FULL_FILE}"; then

        tabix -p vcf "${AF_FILE}"
        AF_SIZE=$(ls -lh "${AF_FILE}" | awk '{print $5}')
        echo "  Full: ${FULL_SIZE} -> AF-only: ${AF_SIZE}"

        # Delete full file to free space for next chromosome
        rm -f "${FULL_FILE}"
        echo "  Deleted full sites VCF"
    else
        echo "  FAILED bcftools annotate on chr${chr}"
        FAILED=$((FAILED + 1))
        rm -f "${AF_FILE}" "${FULL_FILE}"
    fi
    echo ""
done

# Clean up temp dir
rmdir "${TMPDIR_DL}" 2>/dev/null || true

echo "=== gnomAD AF-Only Extraction Complete ==="
echo "End: $(date)"
echo ""

if [[ ${FAILED} -gt 0 ]]; then
    echo "WARNING: ${FAILED} chromosome(s) failed. Re-run to retry."
fi

echo "AF-only files:"
ls -lh "${OUTDIR}"/gnomad.v4.1.af_only.chr*.vcf.bgz 2>/dev/null
echo ""
TOTAL_SIZE=$(du -sh "${OUTDIR}" | awk '{print $1}')
echo "Total gnomAD AF-only size: ${TOTAL_SIZE}"
echo ""
echo "Fields available in AF-only VCFs:"
echo "  AF, AF_grpmax, AC, AN, nhomalt, fafmax_faf95_max"
echo ""
echo "Next: Update run_triage.sh gnomAD paths and field names"
