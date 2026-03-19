#!/usr/bin/env bash
set -euo pipefail

# run_igv_reports.sh - Generate igv-reports HTML for variant validation
#
# Creates an interactive HTML report with embedded BAM reads for visual
# inspection of tiered variants. Supports read strand coloring, MAPQ display,
# and variant allele fraction.
#
# Prerequisites:
#   pip install igv-reports
#
# Usage:
#   ./run_igv_reports.sh -v <tiered_variants.tsv> -b <minibam.bam> [-g <genome>] [-o <output.html>]

GENOME="hg19"
OUTPUT="validation_report.html"
BAM=""
VARIANTS=""
FLANKING=200
FASTA=""

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Generate igv-reports HTML for variant validation.

Options:
  -v TSV       Tiered variants TSV file (required)
               Must have columns: CHROM, POS, REF, ALT (minimum)
  -b BAM       Mini-BAM file (required)
  -g GENOME    Genome build: hg19 or hg38 (default: hg19)
  -f FASTA     Reference FASTA (optional, for read mismatches)
  -w INT       Flanking window in bp (default: 200)
  -o FILE      Output HTML file (default: validation_report.html)
  -h           Show this help

Examples:
  # Basic usage with hg19
  $(basename "$0") -v tiered_variants.tsv -b minibam.bam

  # With reference FASTA for mismatch display
  $(basename "$0") -v tiered_variants.tsv -b minibam.bam -f human_g1k_v37.fasta

  # GRCh38 build
  $(basename "$0") -v tiered_variants.tsv -b minibam.bam -g hg38
EOF
    exit 1
}

while getopts "v:b:g:f:w:o:h" opt; do
    case "$opt" in
        v) VARIANTS="$OPTARG" ;;
        b) BAM="$OPTARG" ;;
        g) GENOME="$OPTARG" ;;
        f) FASTA="$OPTARG" ;;
        w) FLANKING="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "$VARIANTS" || -z "$BAM" ]]; then
    echo "ERROR: both -v and -b are required" >&2
    usage
fi

# Check igv-reports is installed
if ! command -v create_report &>/dev/null; then
    echo "ERROR: igv-reports not installed." >&2
    echo "Install with: pip install igv-reports" >&2
    exit 1
fi

if [[ ! -f "$VARIANTS" ]]; then
    echo "ERROR: Variants file not found: $VARIANTS" >&2
    exit 1
fi

if [[ ! -f "$BAM" ]]; then
    echo "ERROR: BAM file not found: $BAM" >&2
    exit 1
fi

# Check BAM index
if [[ ! -f "${BAM}.bai" && ! -f "${BAM%.bam}.bai" ]]; then
    echo "Indexing BAM..."
    samtools index "$BAM"
fi

NUM_VARIANTS=$(tail -n +2 "$VARIANTS" | wc -l | xargs)
echo ""
echo "=== IGV Report Generation ==="
echo "Variants:  $VARIANTS ($NUM_VARIANTS variants)"
echo "BAM:       $BAM"
echo "Genome:    $GENOME"
echo "Flanking:  ${FLANKING}bp"
echo "Output:    $OUTPUT"
echo ""

# Convert TSV to VCF-like format for igv-reports if needed
# igv-reports can handle various formats; we'll convert to a BED-like table

# Detect if input is VCF or TSV
FIRST_LINE=$(head -1 "$VARIANTS")
if [[ "$FIRST_LINE" == "##"* || "$FIRST_LINE" == "#CHROM"* ]]; then
    INPUT_FORMAT="vcf"
else
    INPUT_FORMAT="tsv"
fi

echo "Detected input format: $INPUT_FORMAT"

# Build igv-reports command
CMD="create_report"
CMD+=" $VARIANTS"
CMD+=" --genome $GENOME"
CMD+=" --flanking $FLANKING"
CMD+=" --tracks $BAM"
CMD+=" --output $OUTPUT"

if [[ -n "$FASTA" ]]; then
    if [[ ! -f "$FASTA" ]]; then
        echo "ERROR: FASTA file not found: $FASTA" >&2
        exit 1
    fi
    CMD+=" --fasta $FASTA"
fi

# For TSV input, specify relevant columns
if [[ "$INPUT_FORMAT" == "tsv" ]]; then
    # Detect column names for sequence (CHROM), begin (POS), end (POS+len(REF))
    # igv-reports needs --sequence, --begin, --end for non-VCF
    HEADER="$FIRST_LINE"

    # Find column names (case-insensitive)
    SEQ_COL=$(echo "$HEADER" | tr '\t' '\n' | grep -i "^chrom$\|^#chrom$\|^chr$\|^chromosome$" | head -1 || true)
    BEGIN_COL=$(echo "$HEADER" | tr '\t' '\n' | grep -i "^pos$\|^position$\|^start$" | head -1 || true)
    END_COL=$(echo "$HEADER" | tr '\t' '\n' | grep -i "^end$\|^stop$" | head -1 || true)

    if [[ -n "$SEQ_COL" ]]; then
        CMD+=" --sequence $SEQ_COL"
    fi
    if [[ -n "$BEGIN_COL" ]]; then
        CMD+=" --begin $BEGIN_COL"
    fi
    if [[ -n "$END_COL" ]]; then
        CMD+=" --end $END_COL"
    elif [[ -n "$BEGIN_COL" ]]; then
        # If no END column, use BEGIN (igv-reports handles single-position)
        CMD+=" --end $BEGIN_COL"
    fi
fi

echo "Running: $CMD"
echo ""

eval "$CMD"

if [[ -f "$OUTPUT" ]]; then
    SIZE=$(ls -lh "$OUTPUT" | awk '{print $5}')
    echo ""
    echo "=== Done ==="
    echo "Report:    $OUTPUT ($SIZE)"
    echo "Variants:  $NUM_VARIANTS"
    echo ""
    echo "Validation criteria for each variant:"
    echo "  1. Bidirectional strand support (reads on both + and - strands)"
    echo "  2. No flanking mismatch clusters (artifact signature)"
    echo "  3. At least 3 alt-supporting reads"
    echo "  4. MAPQ > 0 for supporting reads"
    echo ""
    echo "Open in browser: open $OUTPUT"
else
    echo "ERROR: Report generation failed" >&2
    exit 1
fi
