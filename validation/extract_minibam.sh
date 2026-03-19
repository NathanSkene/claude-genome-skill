#!/usr/bin/env bash
set -euo pipefail

# extract_minibam.sh - Extract mini-BAM regions around variant positions
#
# Extracts reads from a full BAM file for a set of variant positions,
# with configurable flanking region. Designed for IGV visual validation
# of tiered variants from WGS re-analysis.
#
# Usage:
#   ./extract_minibam.sh -b <full.bam> -v <variants.bed> [-f 500] [-o output_dir]
#   ./extract_minibam.sh -b <full.bam> -t <tiered_variants.tsv> [-f 500] [-o output_dir]
#
# After extraction, transfer to local Mac:
#   scp user@hpc:output_dir/minibam.bam* /local/path/

FLANK=500
OUTDIR="./minibam_output"
BAM=""
BED=""
TSV=""

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Extract mini-BAM regions from a full BAM for IGV validation.

Options:
  -b BAM       Input BAM file (required)
  -v BED       BED file of variant positions (chr start end)
  -t TSV       Tiered variants TSV (auto-generates BED from CHROM/POS columns)
  -f INT       Flanking bases each side (default: 500)
  -o DIR       Output directory (default: ./minibam_output)
  -h           Show this help

One of -v or -t is required.

Examples:
  # From BED file
  $(basename "$0") -b aligned.bam -v positions.bed -f 500

  # From tiered variants TSV (expects CHROM and POS columns)
  $(basename "$0") -b aligned.bam -t tiered_variants.tsv

  # Transfer to Mac afterwards
  scp hpc:minibam_output/minibam.bam* ~/Projects/GenomeAnalysis/validation/
EOF
    exit 1
}

while getopts "b:v:t:f:o:h" opt; do
    case "$opt" in
        b) BAM="$OPTARG" ;;
        v) BED="$OPTARG" ;;
        t) TSV="$OPTARG" ;;
        f) FLANK="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "$BAM" ]]; then
    echo "ERROR: -b BAM file is required" >&2
    usage
fi

if [[ -z "$BED" && -z "$TSV" ]]; then
    echo "ERROR: one of -v (BED) or -t (TSV) is required" >&2
    usage
fi

# Check dependencies
for cmd in samtools; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "ERROR: $cmd not found in PATH" >&2
        exit 1
    fi
done

# Check BAM exists and is indexed
if [[ ! -f "$BAM" ]]; then
    echo "ERROR: BAM file not found: $BAM" >&2
    exit 1
fi

if [[ ! -f "${BAM}.bai" && ! -f "${BAM%.bam}.bai" ]]; then
    echo "WARNING: BAM index not found, creating..."
    samtools index "$BAM"
fi

mkdir -p "$OUTDIR"

# Generate BED from TSV if needed
if [[ -n "$TSV" ]]; then
    if [[ ! -f "$TSV" ]]; then
        echo "ERROR: TSV file not found: $TSV" >&2
        exit 1
    fi

    BED="${OUTDIR}/variant_positions.bed"
    echo "Generating BED from TSV: $TSV"

    # Auto-detect CHROM and POS column indices from header
    HEADER=$(head -1 "$TSV")
    CHROM_COL=$(echo "$HEADER" | tr '\t' '\n' | grep -n -i "^chrom$\|^#chrom$\|^chr$" | head -1 | cut -d: -f1)
    POS_COL=$(echo "$HEADER" | tr '\t' '\n' | grep -n -i "^pos$\|^position$\|^start$" | head -1 | cut -d: -f1)

    if [[ -z "$CHROM_COL" || -z "$POS_COL" ]]; then
        echo "ERROR: Could not find CHROM and POS columns in TSV header" >&2
        echo "Header: $HEADER" >&2
        exit 1
    fi

    # Generate BED with flanking regions
    tail -n +2 "$TSV" | awk -v cc="$CHROM_COL" -v pc="$POS_COL" -v flank="$FLANK" \
        'BEGIN{OFS="\t"} {
            chrom = $cc
            pos = $pc
            start = pos - flank
            if (start < 0) start = 0
            end = pos + flank
            print chrom, start, end
        }' | sort -k1,1 -k2,2n | bedtools merge -i - > "$BED" 2>/dev/null || \
        sort -k1,1 -k2,2n > "$BED"

    echo "Generated BED with $(wc -l < "$BED") merged regions"
fi

if [[ ! -f "$BED" ]]; then
    echo "ERROR: BED file not found: $BED" >&2
    exit 1
fi

NUM_REGIONS=$(wc -l < "$BED")
echo ""
echo "=== Mini-BAM Extraction ==="
echo "BAM:      $BAM"
echo "Regions:  $NUM_REGIONS"
echo "Flank:    ${FLANK}bp"
echo ""

# Estimate output size
# Average WGS coverage ~30x, ~150bp reads, ~2 bytes/base in BAM
# Region size * coverage * read_len_factor * BAM_overhead
TOTAL_BP=$(awk '{sum += ($3 - $2)} END {print sum}' "$BED")
EST_SIZE_MB=$(echo "$TOTAL_BP" | awk '{printf "%.1f", ($1 * 30 * 4) / 1048576}')
echo "Total region size: ${TOTAL_BP} bp"
echo "Estimated output:  ~${EST_SIZE_MB} MB (at ~30x coverage)"
echo ""

read -r -p "Proceed with extraction? [Y/n] " CONFIRM
CONFIRM="${CONFIRM:-Y}"
if [[ ! "$CONFIRM" =~ ^[Yy]$ ]]; then
    echo "Aborted."
    exit 0
fi

OUTBAM="${OUTDIR}/minibam.bam"

echo "Extracting reads..."
samtools view -b -h -L "$BED" -o "$OUTBAM" "$BAM"

echo "Sorting..."
samtools sort -o "${OUTBAM%.bam}.sorted.bam" "$OUTBAM"
mv "${OUTBAM%.bam}.sorted.bam" "$OUTBAM"

echo "Indexing..."
samtools index "$OUTBAM"

# Report actual size
ACTUAL_SIZE=$(ls -lh "$OUTBAM" | awk '{print $5}')
echo ""
echo "=== Done ==="
echo "Mini-BAM:  $OUTBAM ($ACTUAL_SIZE)"
echo "Index:     ${OUTBAM}.bai"
echo "Regions:   $NUM_REGIONS"
echo ""
echo "Transfer to local Mac:"
echo "  scp $(hostname):$(realpath "$OUTBAM")* ~/Projects/GenomeAnalysis/validation/"
