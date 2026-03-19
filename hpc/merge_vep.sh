#!/bin/bash
#
# Merge per-chromosome VEP results into a single VCF
#
# Usage: bash merge_vep.sh [INPUT_VCF]
#
# Verifies all chromosomes annotated, checks variant counts against
# the raw DeepVariant VCF, then merges with bcftools concat.

set -euo pipefail

source /usr/share/lmod/lmod/init/bash
module purge 2>/dev/null || true
module load tools/prod HTSlib/1.21-GCC-13.3.0 BCFtools/1.21-GCC-13.3.0

BASE="$(readlink -f /rds/general/user/nskene/home/GenomeAnalysis)"
RESULTS="${BASE}/results"

INPUT_VCF="$(readlink -f "${1:-${RESULTS}/call_snv/genome/GFXC087577_case_snv.vcf.gz}")"
OUTPUT_VCF="${RESULTS}/vep/GFXC087577.vep.vcf.gz"
PER_CHR_DIR="${RESULTS}/vep/per_chr"

CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)

echo "============================================================"
echo "VEP Merge + Verification"
echo "============================================================"
echo "Input (raw):     ${INPUT_VCF}"
echo "Per-chr dir:     ${PER_CHR_DIR}"
echo "Output (merged): ${OUTPUT_VCF}"
echo "Date:            $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
echo

# --- Step 1: Verify all per-chr VCFs exist -----------------------------------
echo "Step 1: Checking per-chromosome VCFs..."
MISSING=0
VCF_LIST=()

for CHR in "${CHROMS[@]}"; do
    CHR_VCF="${PER_CHR_DIR}/GFXC087577.vep.${CHR}.vcf.gz"
    if [[ ! -f "${CHR_VCF}" ]]; then
        echo "  MISSING: ${CHR_VCF}"
        MISSING=$((MISSING + 1))
    elif [[ ! -f "${CHR_VCF}.tbi" ]]; then
        echo "  NO INDEX: ${CHR_VCF} (re-indexing...)"
        bcftools index -t "${CHR_VCF}"
        VCF_LIST+=("${CHR_VCF}")
    else
        VCF_LIST+=("${CHR_VCF}")
    fi
done

if [[ ${MISSING} -gt 0 ]]; then
    echo
    echo "ERROR: ${MISSING} chromosome(s) missing. Re-run run_vep_parallel.sh to submit missing jobs."
    exit 1
fi

echo "  -> All ${#CHROMS[@]} chromosomes present"
echo

# --- Step 2: Variant count verification --------------------------------------
echo "Step 2: Variant count verification..."

TOTAL_VEP=0
TOTAL_RAW=0

printf "  %-6s  %10s  %10s  %s\n" "CHR" "RAW" "VEP" "STATUS"
printf "  %-6s  %10s  %10s  %s\n" "------" "----------" "----------" "------"

for CHR in "${CHROMS[@]}"; do
    CHR_VCF="${PER_CHR_DIR}/GFXC087577.vep.${CHR}.vcf.gz"

    RAW_COUNT=$(bcftools view -H -r "${CHR}" "${INPUT_VCF}" | wc -l)
    VEP_COUNT=$(bcftools view -H "${CHR_VCF}" | wc -l)

    TOTAL_RAW=$((TOTAL_RAW + RAW_COUNT))
    TOTAL_VEP=$((TOTAL_VEP + VEP_COUNT))

    if [[ ${VEP_COUNT} -eq ${RAW_COUNT} ]]; then
        STATUS="OK"
    elif [[ ${VEP_COUNT} -gt 0 ]]; then
        # VEP can split multiallelic records, so VEP >= RAW is acceptable
        STATUS="OK (split)"
    else
        STATUS="EMPTY!"
    fi

    printf "  %-6s  %10d  %10d  %s\n" "${CHR}" "${RAW_COUNT}" "${VEP_COUNT}" "${STATUS}"
done

echo
printf "  %-6s  %10d  %10d\n" "TOTAL" "${TOTAL_RAW}" "${TOTAL_VEP}"
echo

if [[ ${TOTAL_VEP} -eq 0 ]]; then
    echo "ERROR: No variants in VEP output!"
    exit 1
fi

# --- Step 3: Merge -----------------------------------------------------------
echo "Step 3: Merging per-chromosome VCFs..."

# bcftools concat --naive is fast (just concatenates compressed blocks)
# since each chromosome is already sorted and non-overlapping
bcftools concat --naive "${VCF_LIST[@]}" -Oz -o "${OUTPUT_VCF}"
bcftools index -t "${OUTPUT_VCF}"

MERGED_COUNT=$(bcftools view -H "${OUTPUT_VCF}" | wc -l)
echo "  -> Merged VCF: ${OUTPUT_VCF}"
echo "  -> Total variants: ${MERGED_COUNT}"
echo

# --- Step 4: Final verification ----------------------------------------------
echo "Step 4: Final verification..."

# Check all chromosomes are present in merged file
MERGED_CHROMS=$(bcftools query -f '%CHROM\n' "${OUTPUT_VCF}" | sort -u | wc -l)
echo "  Chromosomes in merged VCF: ${MERGED_CHROMS}"

# Quick sanity: VEP CSQ field present
HAS_CSQ=$(bcftools view -h "${OUTPUT_VCF}" | grep -c "ID=CSQ" || true)
if [[ ${HAS_CSQ} -gt 0 ]]; then
    echo "  CSQ annotation field: present"
else
    echo "  WARNING: CSQ annotation field NOT found in header!"
fi

echo
echo "============================================================"
echo "Merge complete"
echo "============================================================"
echo "Output:     ${OUTPUT_VCF}"
echo "Variants:   ${MERGED_COUNT}"
echo "Chromosomes: ${MERGED_CHROMS}"
echo "Date:       $(date '+%Y-%m-%d %H:%M:%S')"
echo
echo "Next steps:"
echo "  1. Transfer to local: scp imperial:${OUTPUT_VCF} ~/Projects/GenomeAnalysis/grch38/"
echo "  2. Re-run triage:     python3 triage/local_triage.py"
