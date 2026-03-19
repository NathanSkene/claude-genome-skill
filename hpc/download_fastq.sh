#!/usr/bin/env bash
set -euo pipefail

# Download BAM from Dante Labs portal for sample ${SAMPLE_ID:-YOUR_SAMPLE}
# Then extract FASTQs using samtools collate/fastq
#
# The Dante Omics portal (genome.danteomics.com) provides presigned
# Cloudflare R2 URLs that are valid for ~1 hour. You must:
#   1. Log in to genome.danteomics.com
#   2. Navigate to profile YOUR_PROFILE_UUID
#   3. Click "Download BAM File" (or FASTQ if available)
#   4. Intercept the presigned URL from the download
#   5. Paste it below or pass as argument

BASEDIR="/rds/general/user/${HPC_USER:-YOUR_USERNAME}/home/GenomeAnalysis"
DLDIR="${BASEDIR}/downloads"
FQDIR="${BASEDIR}/fastq"

mkdir -p "${DLDIR}" "${FQDIR}"

# ──────────────────────────────────────────────────────────────
# Option 1: Download BAM (always available)
# ──────────────────────────────────────────────────────────────
BAM_URL="${1:-}"

if [[ -z "${BAM_URL}" ]]; then
    echo "Usage: $0 <presigned-BAM-URL>"
    echo ""
    echo "Get the URL from genome.danteomics.com:"
    echo "  1. Log in and go to your profile"
    echo "  2. Click 'Download BAM File'"
    echo "  3. Copy the presigned URL (valid ~1 hour)"
    echo ""
    echo "Example:"
    echo "  nohup bash $0 'https://....r2.cloudflarestorage.com/...' > ${DLDIR}/wget.log 2>&1 &"
    exit 1
fi

BAM="${DLDIR}/${SAMPLE_ID:-YOUR_SAMPLE}.bam"

if [[ -f "${BAM}" ]]; then
    BAM_SIZE=$(stat -c%s "${BAM}" 2>/dev/null || stat -f%z "${BAM}")
    echo "BAM already exists: ${BAM} ($(numfmt --to=iec ${BAM_SIZE} 2>/dev/null || echo "${BAM_SIZE} bytes"))"
    echo "Delete it first if you want to re-download."
    exit 1
fi

echo "=== Downloading BAM for ${SAMPLE_ID:-YOUR_SAMPLE} ==="
echo "Output: ${BAM}"
echo "Start: $(date)"
echo ""

wget -O "${BAM}" "${BAM_URL}"

echo ""
echo "Download complete: $(date)"
ls -lh "${BAM}"

echo ""
echo "=== Next Steps ==="
echo "Submit BAM-to-FASTQ conversion job:"
echo "  qsub hpc/bam2fastq.sh"
echo ""
echo "Or if FASTQs are available directly from Dante, download those instead"
echo "and skip the conversion step."
