#!/bin/bash
#PBS -N vep_annotate
#PBS -l select=1:ncpus=8:mem=32gb
#PBS -l walltime=8:00:00
#PBS -q v1_small72
#PBS -o /rds/general/project/genome_analysis/live/GenomeAnalysis/results/vep.out
#PBS -e /rds/general/project/genome_analysis/live/GenomeAnalysis/results/vep.err
#
# Standalone VEP annotation for nf-core/raredisease output
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
INPUT_VCF="$(readlink -f "${1:-${RESULTS}/call_snv/genome/GFXC087577_case_snv.vcf.gz}")"

if [[ ! -f "${INPUT_VCF}" ]]; then
    echo "ERROR: Cannot find DeepVariant VCF. Provide path as argument."
    echo "Usage: qsub -v INPUT_VCF=/path/to/vcf run_vep.sh"
    exit 1
fi

OUTPUT_VCF="${RESULTS}/vep/GFXC087577.vep.vcf.gz"
mkdir -p "${RESULTS}/vep"

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
echo "VEP Annotation"
echo "============================================================"
echo "Input:  ${INPUT_VCF}"
echo "Output: ${OUTPUT_VCF}"
echo "Cache:  ${VEP_CACHE}"
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

# --- Run VEP -----------------------------------------------------------------
echo "Running VEP..."

apptainer exec \
    --bind "${REFS}:${REFS}" \
    --bind "${RESULTS}:${RESULTS}" \
    "${VEP_SIF}" \
    vep \
    --input_file "${INPUT_VCF}" \
    --output_file STDOUT \
    --format vcf \
    --vcf \
    --force_overwrite \
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

echo
echo "VEP annotation complete."
echo "Output: ${OUTPUT_VCF}"
echo "Variants: $(bcftools view -H "${OUTPUT_VCF}" | wc -l)"
echo "End: $(date '+%Y-%m-%d %H:%M:%S')"
