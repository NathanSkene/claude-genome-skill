#!/usr/bin/env bash
set -euo pipefail

# Stage all reference data for nf-core/raredisease on Imperial HPC
# Target: GRCh38 / hg38
# Storage: home RDS (backed up, persistent)

# Clean module environment (avoid .bashrc conflicts with old modules)
source /usr/share/lmod/lmod/init/bash
module purge 2>/dev/null || true
module load tools/prod SAMtools/1.21-GCC-13.3.0

REFDIR="/rds/general/user/${HPC_USER:-YOUR_USERNAME}/home/GenomeAnalysis/references"
mkdir -p "${REFDIR}"

echo "=== Reference Data Staging for nf-core/raredisease ==="
echo "Target directory: ${REFDIR}"
echo ""

# ─────────────────────────────────────────────────────────────
# 1. GRCh38 Reference Genome (GATK bundle / hs38DH)
# ─────────────────────────────────────────────────────────────
GENOME_DIR="${REFDIR}/genome"
mkdir -p "${GENOME_DIR}"

if [[ ! -f "${GENOME_DIR}/Homo_sapiens_assembly38.fasta" ]]; then
    echo "--- Downloading GRCh38 reference genome (GATK bundle) ---"
    cd "${GENOME_DIR}"

    wget -c "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
    wget -c "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
    wget -c "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
    wget -c "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt"

    # BWA-MEM2 index (pre-built saves hours)
    wget -c "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.amb"
    wget -c "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.ann"
    wget -c "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.bwt"
    wget -c "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.pac"
    wget -c "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.sa"

    echo "GRCh38 reference genome: DONE"
else
    echo "GRCh38 reference genome: already staged"
fi

# ─────────────────────────────────────────────────────────────
# 2. VEP Cache (v112, GRCh38)
# ─────────────────────────────────────────────────────────────
VEP_DIR="${REFDIR}/vep"
mkdir -p "${VEP_DIR}"

if [[ ! -d "${VEP_DIR}/homo_sapiens" ]]; then
    echo "--- Downloading VEP cache v112 for GRCh38 (~15 GB) ---"
    cd "${VEP_DIR}"

    wget -c "https://ftp.ensembl.org/pub/release-112/variation/vep/homo_sapiens_vep_112_GRCh38.tar.gz"
    tar xzf homo_sapiens_vep_112_GRCh38.tar.gz
    rm -f homo_sapiens_vep_112_GRCh38.tar.gz

    echo "VEP cache: DONE"
else
    echo "VEP cache: already staged"
fi

# ─────────────────────────────────────────────────────────────
# 3. CADD v1.7 Precomputed Scores (GRCh38)
# ─────────────────────────────────────────────────────────────
CADD_DIR="${REFDIR}/cadd"
mkdir -p "${CADD_DIR}"

if [[ ! -f "${CADD_DIR}/whole_genome_SNVs.tsv.gz" ]]; then
    echo "--- Downloading CADD v1.7 whole-genome SNVs (~80 GB) ---"
    cd "${CADD_DIR}"

    wget -c "https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz"
    wget -c "https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz.tbi"

    echo "CADD SNVs: DONE"
else
    echo "CADD SNVs: already staged"
fi

if [[ ! -f "${CADD_DIR}/gnomad.genomes.r4.0.indel.tsv.gz" ]]; then
    echo "--- Downloading CADD v1.7 indels (~2 GB) ---"
    cd "${CADD_DIR}"

    wget -c "https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz"
    wget -c "https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz.tbi"

    echo "CADD indels: DONE"
else
    echo "CADD indels: already staged"
fi

# ─────────────────────────────────────────────────────────────
# 4. SpliceAI Precomputed Scores (GRCh38)
# ─────────────────────────────────────────────────────────────
SPLICEAI_DIR="${REFDIR}/spliceai"
mkdir -p "${SPLICEAI_DIR}"

if [[ ! -f "${SPLICEAI_DIR}/spliceai_scores.raw.snv.hg38.vcf.gz" ]]; then
    echo "--- Downloading SpliceAI precomputed scores ---"
    echo "NOTE: SpliceAI scores require Illumina BaseSpace download."
    echo "Manual steps:"
    echo "  1. Go to https://basespace.illumina.com/s/otSPW8hnhaZR"
    echo "  2. Download spliceai_scores.raw.snv.hg38.vcf.gz (~17 GB)"
    echo "  3. Download spliceai_scores.raw.snv.hg38.vcf.gz.tbi"
    echo "  4. Download spliceai_scores.raw.indel.hg38.vcf.gz (~4 GB)"
    echo "  5. Download spliceai_scores.raw.indel.hg38.vcf.gz.tbi"
    echo "  6. Place files in: ${SPLICEAI_DIR}/"
    echo ""
    echo "SpliceAI: MANUAL DOWNLOAD REQUIRED"
else
    echo "SpliceAI SNVs: already staged"
fi

# ─────────────────────────────────────────────────────────────
# 5. gnomAD v4 Sites VCF (GRCh38)
# ─────────────────────────────────────────────────────────────
GNOMAD_DIR="${REFDIR}/gnomad"
mkdir -p "${GNOMAD_DIR}"

if [[ ! -f "${GNOMAD_DIR}/gnomad.genomes.v4.1.sites.af-only.vcf.bgz" ]]; then
    echo "--- Downloading gnomAD v4.1 allele frequency VCF ---"
    echo "NOTE: Full gnomAD sites is very large. Downloading AF-only version."
    cd "${GNOMAD_DIR}"

    # AF-only sites VCF (much smaller than full sites)
    wget -c "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr1.vcf.bgz" \
        -O gnomad_chr1_example.vcf.bgz || true

    echo ""
    echo "NOTE: gnomAD v4 full sites VCF is ~700 GB across all chromosomes."
    echo "For nf-core/raredisease, you typically need the AF-only version."
    echo "Consider downloading per-chromosome as needed, or use the bundled"
    echo "gnomAD data from igenomes if available on your HPC."
    echo ""
    echo "Alternative: check if Imperial has gnomAD staged at:"
    echo "  /rds/general/project/reference-data/ (ask HPC team)"
    echo ""
    echo "gnomAD: PARTIAL - check HPC shared references"
else
    echo "gnomAD: already staged"
fi

# ─────────────────────────────────────────────────────────────
# 6. LOFTEE Data (GRCh38)
# ─────────────────────────────────────────────────────────────
LOFTEE_DIR="${REFDIR}/loftee"
mkdir -p "${LOFTEE_DIR}"

if [[ ! -f "${LOFTEE_DIR}/gerp_conservation_scores.homo_sapiens.GRCh38.bw" ]]; then
    echo "--- Downloading LOFTEE data files for GRCh38 ---"
    cd "${LOFTEE_DIR}"

    # GERP conservation scores
    wget -c "https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw"

    # Human ancestor sequence
    wget -c "https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz"
    wget -c "https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.fai"
    wget -c "https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.gzi"

    # PhyloCSF data
    wget -c "https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/loftee.v1.0.4.GRCh38.sql.gz"
    gunzip -k loftee.v1.0.4.GRCh38.sql.gz 2>/dev/null || true

    echo "LOFTEE data: DONE"
else
    echo "LOFTEE data: already staged"
fi

# ─────────────────────────────────────────────────────────────
# 7. pext (GTEx Proportion Expressed Across Transcripts)
# ─────────────────────────────────────────────────────────────
PEXT_DIR="${REFDIR}/pext"
mkdir -p "${PEXT_DIR}"

echo "pext data: SKIPPED (only available as Hail Table, not needed for initial pipeline run)"
echo "  If needed later, export via Hail from: gs://gcp-public-data--gnomad/papers/2019-tx-annotation/pre_computed/"

# ─────────────────────────────────────────────────────────────
# 8. Segmental Duplication BED tracks (GRCh38)
# ─────────────────────────────────────────────────────────────
SEGDUP_DIR="${REFDIR}/segdup"
mkdir -p "${SEGDUP_DIR}"

if [[ ! -f "${SEGDUP_DIR}/GRCh38GenomicSuperDup.bed.gz" ]]; then
    echo "--- Downloading segmental duplication tracks ---"
    cd "${SEGDUP_DIR}"

    wget -c "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz"

    # Convert UCSC table to BED format
    zcat genomicSuperDups.txt.gz | \
        awk -v OFS='\t' '{print $2, $3, $4, $5, $27}' | \
        sort -k1,1 -k2,2n | \
        bgzip > GRCh38GenomicSuperDup.bed.gz

    tabix -p bed GRCh38GenomicSuperDup.bed.gz
    rm -f genomicSuperDups.txt.gz

    echo "Segmental duplications: DONE"
else
    echo "Segmental duplications: already staged"
fi

# ─────────────────────────────────────────────────────────────
# 9. Low-Complexity Region BED tracks (GRCh38)
#    RepeatMasker, TRF, Simple Repeats from UCSC
# ─────────────────────────────────────────────────────────────
LCR_DIR="${REFDIR}/lcr"
mkdir -p "${LCR_DIR}"

if [[ ! -f "${LCR_DIR}/GRCh38_low_complexity_regions.bed.gz" ]]; then
    echo "--- Downloading low-complexity region tracks ---"
    cd "${LCR_DIR}"

    # RepeatMasker
    echo "  Downloading RepeatMasker annotations..."
    wget -c "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz"

    # Simple Repeats (TRF-based)
    echo "  Downloading Simple Repeats (TRF)..."
    wget -c "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz"

    # Convert RepeatMasker to BED
    echo "  Converting RepeatMasker to BED..."
    zcat rmsk.txt.gz | \
        awk -v OFS='\t' '{print $6, $7, $8, $11 "/" $12 "/" $13}' | \
        sort -k1,1 -k2,2n | \
        bgzip > GRCh38_repeatmasker.bed.gz
    tabix -p bed GRCh38_repeatmasker.bed.gz

    # Convert Simple Repeats to BED
    echo "  Converting Simple Repeats to BED..."
    zcat simpleRepeat.txt.gz | \
        awk -v OFS='\t' '{print $2, $3, $4, "SimpleRepeat"}' | \
        sort -k1,1 -k2,2n | \
        bgzip > GRCh38_simple_repeats.bed.gz
    tabix -p bed GRCh38_simple_repeats.bed.gz

    # Merge into combined low-complexity BED
    echo "  Creating merged low-complexity BED..."
    module load BEDTools/2.31.0-GCC-12.3.0 2>/dev/null || true

    if command -v bedtools &>/dev/null; then
        zcat GRCh38_repeatmasker.bed.gz GRCh38_simple_repeats.bed.gz | \
            cut -f1-3 | \
            sort -k1,1 -k2,2n | \
            bedtools merge -i - | \
            bgzip > GRCh38_low_complexity_regions.bed.gz
        tabix -p bed GRCh38_low_complexity_regions.bed.gz
    else
        echo "  WARNING: bedtools not found. Merge step skipped."
        echo "  Run: module load bedtools && bedtools merge ..."
    fi

    rm -f rmsk.txt.gz simpleRepeat.txt.gz

    echo "Low-complexity regions: DONE"
else
    echo "Low-complexity regions: already staged"
fi

# ─────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────
echo ""
echo "=== Reference Data Summary ==="
echo ""

check_staged() {
    local label="$1"
    local path="$2"
    if [[ -e "${path}" ]]; then
        local size
        size=$(du -sh "${path}" 2>/dev/null | cut -f1)
        printf "  %-35s STAGED (%s)\n" "${label}" "${size}"
    else
        printf "  %-35s MISSING\n" "${label}"
    fi
}

check_staged "GRCh38 reference genome" "${GENOME_DIR}/Homo_sapiens_assembly38.fasta"
check_staged "VEP cache v112" "${VEP_DIR}/homo_sapiens"
check_staged "CADD v1.7 SNVs" "${CADD_DIR}/whole_genome_SNVs.tsv.gz"
check_staged "CADD v1.7 indels" "${CADD_DIR}/gnomad.genomes.r4.0.indel.tsv.gz"
check_staged "SpliceAI SNVs" "${SPLICEAI_DIR}/spliceai_scores.raw.snv.hg38.vcf.gz"
check_staged "SpliceAI indels" "${SPLICEAI_DIR}/spliceai_scores.raw.indel.hg38.vcf.gz"
check_staged "gnomAD v4" "${GNOMAD_DIR}/gnomad.genomes.v4.1.sites.af-only.vcf.bgz"
check_staged "LOFTEE data" "${LOFTEE_DIR}/gerp_conservation_scores.homo_sapiens.GRCh38.bw"
echo "  pext data                           SKIPPED (Hail Table only)"
check_staged "Segmental duplications" "${SEGDUP_DIR}/GRCh38GenomicSuperDup.bed.gz"
check_staged "Low-complexity regions" "${LCR_DIR}/GRCh38_low_complexity_regions.bed.gz"

echo ""
echo "Total reference disk usage:"
du -sh "${REFDIR}"
echo ""
echo "NOTE: Expect ~120-150 GB total for all references."
echo "SpliceAI requires manual download from Illumina BaseSpace."
echo "gnomAD v4 full sites is very large - check if Imperial has shared copies."
