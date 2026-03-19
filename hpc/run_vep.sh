#!/bin/bash
#PBS -l select=1:ncpus=8:mem=32gb
#PBS -l walltime=4:00:00
#PBS -q v1_small72
#
# Standalone VEP annotation for nf-core/raredisease output
# Supports whole-genome or per-chromosome mode via CHR variable.
#
# Usage:
#   Whole genome: qsub run_vep.sh
#   Single chr:   qsub -N vep_chr1 -v CHR=chr1 run_vep.sh
#   Parallel:     bash run_vep_parallel.sh  (submits all chromosomes)
#
# Input: DeepVariant VCF from raredisease pipeline
# Output: VEP-annotated VCF ready for triage engine

set -euo pipefail

# --- Environment setup -------------------------------------------------------
source /usr/share/lmod/lmod/init/bash
module purge 2>/dev/null || true
module load tools/prod HTSlib/1.21-GCC-13.3.0 BCFtools/1.21-GCC-13.3.0

# --- Paths -------------------------------------------------------------------
# Resolve symlinks — Apptainer bind mounts require real paths, not symlinks
BASE="$(readlink -f /rds/general/user/nskene/home/GenomeAnalysis)"
RESULTS="${BASE}/results"
REFS="${BASE}/references"
CONTAINERS="${BASE}/containers"

# Input: DeepVariant VCF from nf-core/raredisease call_snv output
INPUT_VCF="$(readlink -f "${INPUT_VCF:-${RESULTS}/call_snv/genome/GFXC087577_case_snv.vcf.gz}")"

if [[ ! -f "${INPUT_VCF}" ]]; then
    echo "ERROR: Cannot find DeepVariant VCF at: ${INPUT_VCF}"
    echo "Usage: qsub -v INPUT_VCF=/path/to/vcf run_vep.sh"
    exit 1
fi

# Per-chromosome mode: CHR variable (e.g. chr1, chrX)
# If CHR is set, annotate only that chromosome; otherwise annotate everything
CHR="${CHR:-}"

if [[ -n "${CHR}" ]]; then
    OUTPUT_VCF="${RESULTS}/vep/per_chr/GFXC087577.vep.${CHR}.vcf.gz"
    mkdir -p "${RESULTS}/vep/per_chr"
else
    OUTPUT_VCF="${RESULTS}/vep/GFXC087577.vep.vcf.gz"
    mkdir -p "${RESULTS}/vep"
fi

# VEP cache and plugin data
VEP_CACHE="${REFS}/vep"
CADD_SNV="${REFS}/cadd/whole_genome_SNVs.tsv.gz"
CADD_INDEL="${REFS}/cadd/gnomad.genomes.r4.0.indel.tsv.gz"
LOFTEE_DIR="${REFS}/loftee"
SPLICEAI_SNV="${REFS}/spliceai/spliceai_scores.raw.snv.hg38.vcf.gz"
SPLICEAI_INDEL="${REFS}/spliceai/spliceai_scores.raw.indel.hg38.vcf.gz"
ALPHAMISSENSE="${REFS}/alphamissense/AlphaMissense_hg38.tsv.gz"
REVEL="${REFS}/revel/new_tabbed_revel_grch38.tsv.gz"
FASTA="${REFS}/genome/Homo_sapiens_assembly38.fasta"

# VEP container — v114 brings MaveDB, PrimateAI-3D, updated LOEUF v4.1, AllOfUs frequencies
VEP_SIF="${CONTAINERS}/ensembl-vep_114.0.sif"

echo "============================================================"
echo "VEP Annotation${CHR:+ (${CHR})}"
echo "============================================================"
echo "Input:  ${INPUT_VCF}"
echo "Output: ${OUTPUT_VCF}"
echo "Cache:  ${VEP_CACHE}"
echo "Chr:    ${CHR:-all}"
echo "Date:   $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
echo

# --- Pull VEP container if not present ---------------------------------------
if [[ ! -f "${VEP_SIF}" ]]; then
    echo "Pulling VEP v114 container..."
    mkdir -p "${CONTAINERS}"
    apptainer pull "${VEP_SIF}" docker://ensemblorg/ensembl-vep:release_114.0
    echo "  -> Done"
    echo
fi

# --- Extract chromosome if per-chr mode --------------------------------------
if [[ -n "${CHR}" ]]; then
    CHR_VCF="${RESULTS}/vep/per_chr/input.${CHR}.vcf.gz"
    echo "Extracting ${CHR} from input VCF..."
    bcftools view -r "${CHR}" "${INPUT_VCF}" -Oz -o "${CHR_VCF}"
    bcftools index -t "${CHR_VCF}"
    VEP_INPUT="${CHR_VCF}"
    VARIANTS_IN=$(bcftools view -H "${CHR_VCF}" | wc -l)
    echo "  -> ${VARIANTS_IN} variants on ${CHR}"
    echo
else
    VEP_INPUT="${INPUT_VCF}"
fi

# --- Run VEP -----------------------------------------------------------------
echo "Running VEP..."

apptainer exec \
    --bind "${REFS}:${REFS}" \
    --bind "${RESULTS}:${RESULTS}" \
    "${VEP_SIF}" \
    vep \
    --input_file "${VEP_INPUT}" \
    --output_file STDOUT \
    --format vcf \
    --vcf \
    --force_overwrite \
    --check_existing \
    --offline \
    --cache \
    --dir_cache "${VEP_CACHE}" \
    --species homo_sapiens \
    --assembly GRCh38 \
    --fasta "${FASTA}" \
    --fork 8 \
    --buffer_size 5000 \
    --everything \
    --allele_number \
    --no_stats \
    --pick \
    --pick_allele_gene \
    --plugin CADD,"${CADD_SNV}","${CADD_INDEL}" \
    --plugin LoF,loftee_path:/opt/vep/Plugins,human_ancestor_fa:"${LOFTEE_DIR}/human_ancestor.fa.gz",gerp_bigwig:"${LOFTEE_DIR}/gerp_conservation_scores.homo_sapiens.GRCh38.bw" \
    --plugin SpliceAI,snv="${SPLICEAI_SNV}",indel="${SPLICEAI_INDEL}" \
    --plugin AlphaMissense,file="${ALPHAMISSENSE}" \
    --plugin REVEL,"${REVEL}" \
    | bcftools view -Oz -o "${OUTPUT_VCF}"

bcftools index -t "${OUTPUT_VCF}"

# Clean up per-chr input extract
if [[ -n "${CHR}" && -f "${CHR_VCF}" ]]; then
    rm -f "${CHR_VCF}" "${CHR_VCF}.tbi"
fi

echo
echo "VEP annotation complete."
echo "Output: ${OUTPUT_VCF}"
echo "Variants: $(bcftools view -H "${OUTPUT_VCF}" | wc -l)"
echo "End: $(date '+%Y-%m-%d %H:%M:%S')"
