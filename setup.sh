#!/bin/bash
# Bootstrap script for claude-genome-skill
# Prepares your genome VCF for use with the skill

set -euo pipefail

# --- Configuration ---
GENOME_DIR="${GENOME_DIR:-$HOME/genome}"

echo "=== Claude Genome Skill Setup ==="
echo "GENOME_DIR: $GENOME_DIR"
echo ""

# --- Check dependencies ---
echo "Checking dependencies..."

check_dep() {
    if ! command -v "$1" &> /dev/null; then
        echo "  ✗ $1 not found"
        echo "    Install with: $2"
        return 1
    else
        echo "  ✓ $1 ($($1 --version 2>&1 | head -1))"
        return 0
    fi
}

MISSING=0
check_dep bcftools "brew install bcftools (macOS) or apt install bcftools (Ubuntu)" || MISSING=1
check_dep tabix "brew install htslib (macOS) or apt install tabix (Ubuntu)" || MISSING=1
check_dep python3 "brew install python3 (macOS) or apt install python3 (Ubuntu)" || MISSING=1

if [ $MISSING -eq 1 ]; then
    echo ""
    echo "Install missing dependencies and re-run this script."
    exit 1
fi

# Check Python packages
echo ""
echo "Checking Python packages..."
python3 -c "import cyvcf2" 2>/dev/null && echo "  ✓ cyvcf2" || { echo "  ✗ cyvcf2 — pip3 install cyvcf2"; MISSING=1; }
python3 -c "from biothings_client import get_client" 2>/dev/null && echo "  ✓ biothings-client" || { echo "  ✗ biothings-client — pip3 install biothings-client"; MISSING=1; }

if [ $MISSING -eq 1 ]; then
    echo ""
    echo "Install missing Python packages:"
    echo "  pip3 install cyvcf2 biothings-client"
    exit 1
fi

# --- Setup genome directory ---
echo ""
echo "Setting up genome directory..."
mkdir -p "$GENOME_DIR"
mkdir -p "$GENOME_DIR/grch38"

# --- Check for VCF files ---
echo ""
echo "Looking for VCF files..."

VCF_FOUND=0
for vcf in "$GENOME_DIR/grch38/genome.grch38.annotated.vcf.gz" \
           "$GENOME_DIR/grch38/genome.grch38.vcf.gz" \
           "$GENOME_DIR/genome.rsid.vcf.gz" \
           "$GENOME_DIR/genome.vcf.gz"; do
    if [ -f "$vcf" ]; then
        echo "  ✓ Found: $vcf"
        VCF_FOUND=1

        # Check for index
        if [ -f "${vcf}.tbi" ]; then
            echo "    ✓ Index exists"
        elif [ -f "${vcf}.csi" ]; then
            echo "    ✓ CSI index exists"
        else
            echo "    ✗ No index — creating..."
            tabix -p vcf "$vcf"
            echo "    ✓ Index created"
        fi
    fi
done

if [ $VCF_FOUND -eq 0 ]; then
    echo "  ✗ No VCF files found in $GENOME_DIR"
    echo ""
    echo "  Place your genome VCF in one of these locations:"
    echo "    $GENOME_DIR/genome.vcf.gz              (GRCh37)"
    echo "    $GENOME_DIR/genome.rsid.vcf.gz         (GRCh37, rsID-annotated)"
    echo "    $GENOME_DIR/grch38/genome.grch38.vcf.gz (GRCh38)"
    echo ""
    echo "  Then re-run this script to index it."
    exit 0
fi

# --- Check for PRS results ---
echo ""
if [ -f "$GENOME_DIR/prs_results.json" ]; then
    N_SCORES=$(python3 -c "import json; d=json.load(open('$GENOME_DIR/prs_results.json')); print(d.get('n_scores', 'unknown'))" 2>/dev/null || echo "?")
    echo "✓ PRS results found ($N_SCORES scores)"
else
    echo "ℹ No PRS results found (optional — compute with pgsc_calc on HPC)"
fi

# --- Add rsIDs (optional) ---
echo ""
if [ -f "$GENOME_DIR/genome.vcf.gz" ] && [ ! -f "$GENOME_DIR/genome.rsid.vcf.gz" ]; then
    echo "Tip: Your VCF doesn't have rsID annotations. To add them:"
    echo ""
    echo "  # Download dbSNP for your build"
    echo "  # GRCh37:"
    echo "  wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz"
    echo "  wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz.tbi"
    echo ""
    echo "  # Annotate"
    echo "  bcftools annotate -a 00-All.vcf.gz -c ID $GENOME_DIR/genome.vcf.gz -Oz -o $GENOME_DIR/genome.rsid.vcf.gz"
    echo "  tabix -p vcf $GENOME_DIR/genome.rsid.vcf.gz"
fi

# --- Environment variable reminder ---
echo ""
echo "=== Setup Complete ==="
echo ""
echo "Add to your shell profile (~/.zshrc or ~/.bashrc):"
echo ""
echo "  export GENOME_DIR=$GENOME_DIR"
echo ""
echo "Install the skill:"
echo ""
echo "  # If not already in skills directory:"
echo "  ln -s $(pwd) ~/.claude/skills/genome"
echo ""
echo "  # Copy trigger rules for auto-detection:"
echo "  cp genome-triggers.md ~/.claude/rules/"
echo ""
echo "Test it:"
echo ""
echo "  python3 genome_query.py stats"
echo "  python3 genome_query.py variant rs1805007"
