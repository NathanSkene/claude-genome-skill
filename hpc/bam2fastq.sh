#!/usr/bin/env bash
#PBS -N bam2fastq
#PBS -l walltime=06:00:00,select=1:ncpus=8:mem=128gb
#PBS -o /rds/general/user/${HPC_USER:-YOUR_USERNAME}/home/GenomeAnalysis/downloads/bam2fastq.out
#PBS -e /rds/general/user/${HPC_USER:-YOUR_USERNAME}/home/GenomeAnalysis/downloads/bam2fastq.err

set -euo pipefail

source /usr/share/lmod/lmod/init/bash
module purge 2>/dev/null || true
module load tools/prod SAMtools/1.21-GCC-13.3.0

BASEDIR=/rds/general/user/${HPC_USER:-YOUR_USERNAME}/home/GenomeAnalysis
BAM="/rds/general/user/${HPC_USER:-YOUR_USERNAME}/ephemeral/${SAMPLE_ID:-YOUR_SAMPLE}.bam"
R1="${BASEDIR}/fastq/${SAMPLE_ID:-YOUR_SAMPLE}_R1.fastq.gz"
R2="${BASEDIR}/fastq/${SAMPLE_ID:-YOUR_SAMPLE}_R2.fastq.gz"

echo "=== BAM to FASTQ conversion ==="
echo "Start: $(date)"
echo "BAM: ${BAM}"
echo "Output R1: ${R1}"
echo "Output R2: ${R2}"
echo ""

# Verify BAM exists and is not truncated
if [[ ! -f "${BAM}" ]]; then
    echo "ERROR: BAM file not found: ${BAM}"
    exit 1
fi

BAM_SIZE=$(stat -c%s "${BAM}" 2>/dev/null || stat -f%z "${BAM}")
echo "BAM size: ${BAM_SIZE} bytes ($(numfmt --to=iec ${BAM_SIZE} 2>/dev/null || echo "${BAM_SIZE}"))"

if [[ ${BAM_SIZE} -lt 1000000000 ]]; then
    echo "ERROR: BAM seems too small (<1GB). Download may be incomplete."
    exit 1
fi

# collate (name-sort without full sort overhead) then extract paired FASTQs
# -u: uncompressed output from collate (piped, so no disk)
# -O: output to stdout
# -n: don't append /1 /2 to read names
# -0 /dev/null: discard unpaired reads
# -s /dev/null: discard singleton reads
# --threads 7: use 7 additional threads (8 total with main thread)

# Use ephemeral for temp files — NOT PBS $TMPDIR which gets cleaned up
# before the pipe fully flushes, causing truncated output
SCRATCH="/rds/general/user/${HPC_USER:-YOUR_USERNAME}/ephemeral/bam2fastq_tmp"
mkdir -p "${SCRATCH}"
echo "Temp dir: ${SCRATCH}"
echo "Running samtools collate | samtools fastq..."
samtools collate -u -O "${BAM}" "${SCRATCH}/collate_tmp_$$" | \
samtools fastq -1 "${R1}" \
               -2 "${R2}" \
               -0 /dev/null -s /dev/null \
               -n --threads 7

echo ""
echo "=== Verification ==="

# Check files exist and have reasonable size
for f in "${R1}" "${R2}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: Output file missing: ${f}"
        exit 1
    fi
    SIZE=$(stat -c%s "${f}" 2>/dev/null || stat -f%z "${f}")
    echo "$(basename ${f}): ${SIZE} bytes ($(numfmt --to=iec ${SIZE} 2>/dev/null || echo "${SIZE}"))"
done

# Spot-check: first 8 lines of R1 should be 2 valid FASTQ records
echo ""
echo "First 2 records of R1:"
zcat "${R1}" | head -8

# Line count check (should be divisible by 4)
echo ""
echo "Counting lines in R1 (this may take a while)..."
R1_LINES=$(zcat "${R1}" | wc -l)
echo "R1 lines: ${R1_LINES}"
if (( R1_LINES % 4 != 0 )); then
    echo "WARNING: R1 line count not divisible by 4!"
else
    echo "R1 reads: $((R1_LINES / 4))"
fi

# Generate md5sums
echo ""
echo "Generating md5sums..."
md5sum "${R1}" "${R2}" | tee "${BASEDIR}/fastq/md5sums.txt"

echo ""
echo "Done: $(date)"
echo ""
echo "Next steps:"
echo "  1. Verify read counts look reasonable (~900M-1.2B reads per file for 30x WGS)"
echo "  2. Update fastq/README.md with actual sizes and checksums"
echo "  3. Delete BAM to reclaim space: rm ${BAM}"
