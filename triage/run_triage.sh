#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# WGS Variant Triage Pipeline
# Artifact-aware filtering and tiering using vcfanno + slivar
#
# Usage: ./run_triage.sh <annotated_vcf>
#
# Input: VEP-annotated VCF from nf-core/raredisease
# Output: Tiered variant TSVs + artifact reports + module reports
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"

# --- HPC Module Setup --------------------------------------------------------

if [[ -f /usr/share/lmod/lmod/init/bash ]]; then
    source /usr/share/lmod/lmod/init/bash
    module purge 2>/dev/null || true
    module load tools/prod HTSlib/1.21-GCC-13.3.0 BCFtools/1.21-GCC-13.3.0
fi

# --- Configuration -----------------------------------------------------------
# CONFIGURE: Update these paths for your HPC environment

INPUT_VCF="${1:?Usage: $0 <annotated_vcf>}"
OUTPUT_DIR="${SCRIPT_DIR}/output"

# Reference paths (HPC)
REF_DIR="/rds/general/user/${HPC_USER:-YOUR_USERNAME}/home/GenomeAnalysis/references/triage"
GNOMAD_DIR="${REF_DIR}/gnomad"

VCFANNO_TOML="${SCRIPT_DIR}/vcfanno.toml"
SLIVAR_JS="${SCRIPT_DIR}/slivar_filters.js"
KOSMOS_SCRIPT="${SCRIPT_DIR}/kosmos_check.sh"

# Tool paths — prefer Apptainer containers if pulled
CONTAINER_DIR="/rds/general/user/${HPC_USER:-YOUR_USERNAME}/home/GenomeAnalysis/containers"
if [[ -f "${CONTAINER_DIR}/vcfanno_0.3.7.sif" ]]; then
    VCFANNO="apptainer exec ${CONTAINER_DIR}/vcfanno_0.3.7.sif vcfanno"
else
    VCFANNO="${VCFANNO:-vcfanno}"
fi
if [[ -f "${CONTAINER_DIR}/slivar_0.3.2.sif" ]]; then
    SLIVAR="apptainer exec ${CONTAINER_DIR}/slivar_0.3.2.sif slivar"
else
    SLIVAR="${SLIVAR:-slivar}"
fi
BCFTOOLS="${BCFTOOLS:-bcftools}"

# Output files
ANNOTATED_VCF="${OUTPUT_DIR}/context_annotated.vcf.gz"
FLAGGED_VCF="${OUTPUT_DIR}/flagged.vcf.gz"
TIERED_ALL="${OUTPUT_DIR}/tiered_variants.tsv"
TIER1_TSV="${OUTPUT_DIR}/tier1.tsv"
TIER2_TSV="${OUTPUT_DIR}/tier2.tsv"
TIER3_TSV="${OUTPUT_DIR}/tier3.tsv"
CLUSTER_VCF="${OUTPUT_DIR}/cluster_flagged.vcf.gz"
SUMMARY="${OUTPUT_DIR}/triage_summary.txt"

mkdir -p "${OUTPUT_DIR}"

echo "============================================================"
echo "WGS Variant Triage Pipeline"
echo "============================================================"
echo "Input:  ${INPUT_VCF}"
echo "Output: ${OUTPUT_DIR}/"
echo "Date:   $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
echo

# =============================================================================
# Step 0: gnomAD annotation (per-chromosome bcftools annotate)
# =============================================================================

echo "[Step 0/8] Annotating with gnomAD v4.1 allele frequencies..."

GNOMAD_INPUT="${INPUT_VCF}"
GNOMAD_OUTPUT="${OUTPUT_DIR}/gnomad_annotated.vcf.gz"

if [[ -d "${GNOMAD_DIR}" ]]; then
    # Build bcftools annotate command for per-chr gnomAD files
    # We annotate chromosome by chromosome and pipe through
    GNOMAD_TMPVCF="${OUTPUT_DIR}/gnomad_tmp.vcf.gz"
    cp "${GNOMAD_INPUT}" "${GNOMAD_TMPVCF}"
    ${BCFTOOLS} index -t "${GNOMAD_TMPVCF}" 2>/dev/null || true

    ANNOTATED_CHROMS=0
    for chr in $(seq 1 22) X Y; do
        GNOMAD_CHR="${GNOMAD_DIR}/gnomad.v4.1.af_only.chr${chr}.vcf.bgz"
        if [[ -f "${GNOMAD_CHR}" && -f "${GNOMAD_CHR}.tbi" ]]; then
            echo "  chr${chr}..."
            GNOMAD_NEXT="${OUTPUT_DIR}/gnomad_chr${chr}.vcf.gz"
            ${BCFTOOLS} annotate \
                -a "${GNOMAD_CHR}" \
                -c INFO/AF,INFO/AF_grpmax,INFO/nhomalt,INFO/fafmax_faf95_max \
                -r "chr${chr}" \
                -Oz -o "${GNOMAD_NEXT}" \
                "${GNOMAD_TMPVCF}"
            ${BCFTOOLS} index -t "${GNOMAD_NEXT}"
            mv "${GNOMAD_NEXT}" "${GNOMAD_TMPVCF}"
            mv "${GNOMAD_NEXT}.tbi" "${GNOMAD_TMPVCF}.tbi"
            ANNOTATED_CHROMS=$((ANNOTATED_CHROMS + 1))
        fi
    done
    mv "${GNOMAD_TMPVCF}" "${GNOMAD_OUTPUT}"
    mv "${GNOMAD_TMPVCF}.tbi" "${GNOMAD_OUTPUT}.tbi"
    echo "  -> Annotated ${ANNOTATED_CHROMS} chromosomes with gnomAD AF"
    VCFANNO_INPUT="${GNOMAD_OUTPUT}"
else
    echo "  [WARN] gnomAD directory not found: ${GNOMAD_DIR}"
    echo "  -> Proceeding without gnomAD annotation"
    VCFANNO_INPUT="${INPUT_VCF}"
fi
echo

# =============================================================================
# Step 1: vcfanno - Add genomic context annotations
# =============================================================================

echo "[Step 1/8] Running vcfanno to add genomic context annotations..."

${VCFANNO} -p 4 "${VCFANNO_TOML}" "${VCFANNO_INPUT}" \
    | ${BCFTOOLS} view -Oz -o "${ANNOTATED_VCF}"
${BCFTOOLS} index -t "${ANNOTATED_VCF}"

echo "  -> Added: SEGDUP, LCR, TRF, SIMPLE_REPEAT"
echo

# =============================================================================
# Step 2: Positional clustering detection (pre-slivar)
# =============================================================================

echo "[Step 2/8] Detecting positional clustering of indels..."

# Flag indels where another indel exists within 10bp
# Uses bcftools to extract indels and python to detect clusters
${BCFTOOLS} query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' \
    -i 'STRLEN(REF)!=STRLEN(ALT)' "${ANNOTATED_VCF}" \
    | python3 -c "
import sys
from collections import defaultdict

# Read all indels
indels = []
for line in sys.stdin:
    parts = line.strip().split('\t')
    chrom, pos = parts[0], int(parts[1])
    indels.append((chrom, pos, parts[2], parts[3], parts[4]))

# Find clusters: multiple indels within 10bp
clustered = set()
for i, (c1, p1, r1, a1, id1) in enumerate(indels):
    for j, (c2, p2, r2, a2, id2) in enumerate(indels):
        if i != j and c1 == c2 and abs(p1 - p2) <= 10:
            clustered.add((c1, p1))
            clustered.add((c2, p2))

# Output as BED for annotation
for chrom, pos in sorted(clustered):
    print(f'{chrom}\t{pos-1}\t{pos}\tCLUSTER_INDEL')
" > "${OUTPUT_DIR}/clustered_indels.bed"

CLUSTER_COUNT=$(wc -l < "${OUTPUT_DIR}/clustered_indels.bed" | tr -d ' ')
echo "  -> Found ${CLUSTER_COUNT} indels in clusters (within 10bp of another indel)"

# If clusters exist, annotate the VCF with cluster flag
if [ "${CLUSTER_COUNT}" -gt 0 ]; then
    bgzip -c "${OUTPUT_DIR}/clustered_indels.bed" > "${OUTPUT_DIR}/clustered_indels.bed.gz"
    tabix -p bed "${OUTPUT_DIR}/clustered_indels.bed.gz"

    # Add CLUSTER_INDEL flag via bcftools annotate
    echo '##INFO=<ID=CLUSTER_INDEL,Number=0,Type=Flag,Description="Indel within 10bp of another indel (potential artifact)">' \
        > "${OUTPUT_DIR}/cluster_header.txt"

    ${BCFTOOLS} annotate \
        -a "${OUTPUT_DIR}/clustered_indels.bed.gz" \
        -c CHROM,FROM,TO,CLUSTER_INDEL \
        -h "${OUTPUT_DIR}/cluster_header.txt" \
        -Oz -o "${CLUSTER_VCF}" \
        "${ANNOTATED_VCF}"
    ${BCFTOOLS} index -t "${CLUSTER_VCF}"
    WORKING_VCF="${CLUSTER_VCF}"
else
    WORKING_VCF="${ANNOTATED_VCF}"
fi
echo

# =============================================================================
# Step 3: slivar - Apply artifact flags and tiering
# =============================================================================

echo "[Step 3/8] Running slivar for artifact flagging and tiering..."

${SLIVAR} expr \
    --js "${SLIVAR_JS}" \
    --vcf "${WORKING_VCF}" \
    --info "cluster_artifact(variant)" \
    --info "ab_skew(variant)" \
    --info "low_quality(variant)" \
    --info "biological_implausibility(variant)" \
    --info "tier1_diagnostic(variant)" \
    --info "tier2_strong(variant)" \
    --info "tier3_moderate(variant)" \
    --out-vcf "${FLAGGED_VCF}"

${BCFTOOLS} index -t "${FLAGGED_VCF}"

echo "  -> Artifact flags: CLUSTER_ARTIFACT, AB_SKEW, LOW_QUAL, BIOLOGICAL_IMPLAUSIBILITY"
echo "  -> Tier assignments: TIER1_DIAGNOSTIC, TIER2_STRONG, TIER3_MODERATE"
echo

# =============================================================================
# Step 4: Generate tiered TSV output
# =============================================================================

echo "[Step 4/8] Generating tiered variant tables..."

TSV_HEADER="CHROM\tPOS\tREF\tALT\tGENE\tCONSEQUENCE\tTIER\tFLAGS\tCLINVAR\tGNOMAD_AF\tCADD\tSPLICEAI\tLOFTEE\tALPHAMISSENSE\tREVEL\tPEXT\tZYGOSITY\tDP\tGQ\tAB"

# Extract fields from flagged VCF into TSV
# Uses bcftools query with computed fields
${BCFTOOLS} query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CSQ_SYMBOL\t%INFO/CSQ_Consequence\t.\t.\t%INFO/CSQ_ClinVar_CLNSIG\t%INFO/gnomAD_AF\t%INFO/CSQ_CADD_PHRED\t%INFO/CSQ_SpliceAI_DS_max\t%INFO/CSQ_LoF\t%INFO/CSQ_am_class\t%INFO/CSQ_REVEL\t%INFO/CSQ_pext\t[%GT]\t[%DP]\t[%GQ]\t[%AD]\n' \
    "${FLAGGED_VCF}" 2>/dev/null \
    | python3 -c "
import sys

print('CHROM\tPOS\tREF\tALT\tGENE\tCONSEQUENCE\tTIER\tFLAGS\tCLINVAR\tGNOMAD_AF\tCADD\tSPLICEAI\tLOFTEE\tALPHAMISSENSE\tREVEL\tPEXT\tZYGOSITY\tDP\tGQ\tAB')

for line in sys.stdin:
    fields = line.strip().split('\t')
    if len(fields) < 20:
        continue

    chrom, pos, ref, alt = fields[0:4]
    gene = fields[4] if fields[4] != '.' else ''
    consequence = fields[5] if fields[5] != '.' else ''
    clinvar = fields[8] if fields[8] != '.' else ''
    gnomad_af = fields[9] if fields[9] != '.' else ''
    cadd = fields[10] if fields[10] != '.' else ''
    spliceai = fields[11] if fields[11] != '.' else ''
    loftee = fields[12] if fields[12] != '.' else ''
    alphamissense = fields[13] if fields[13] != '.' else ''
    revel = fields[14] if fields[14] != '.' else ''
    pext = fields[15] if fields[15] != '.' else ''
    gt = fields[16]
    dp = fields[17]
    gq = fields[18]
    ad_raw = fields[19]

    # Compute zygosity from GT
    if gt in ('1/1', '1|1'):
        zygosity = 'HOM_ALT'
    elif gt in ('0/1', '0|1', '1|0'):
        zygosity = 'HET'
    elif gt in ('0/0', '0|0'):
        zygosity = 'HOM_REF'
    else:
        zygosity = gt

    # Compute allele balance from AD
    ab = ''
    if ad_raw and ad_raw != '.':
        try:
            ads = [int(x) for x in ad_raw.split(',')]
            if len(ads) >= 2 and sum(ads) > 0:
                ab = f'{ads[1]/sum(ads):.3f}'
        except ValueError:
            pass

    # Assign tier (check slivar INFO flags)
    # Since we extracted from the flagged VCF, we re-derive tier from evidence
    tier = ''
    if clinvar and ('Pathogenic' in clinvar or 'Likely_pathogenic' in clinvar):
        tier = 'TIER1'
    elif loftee == 'HC' or (spliceai and spliceai != '' and float(spliceai) > 0.8 if spliceai else False):
        tier = 'TIER2'
    elif (cadd and cadd != '' and revel and revel != ''):
        try:
            if float(cadd) > 20 and float(revel) > 0.5:
                tier = 'TIER3'
        except ValueError:
            pass

    # Build flags list
    flags = []
    if ad_raw and ab:
        try:
            ab_val = float(ab)
            if zygosity == 'HET' and (ab_val < 0.25 or ab_val > 0.75):
                flags.append('AB_SKEW')
        except ValueError:
            pass
    if dp and dp != '.':
        try:
            if int(dp) < 15:
                flags.append('LOW_DP')
        except ValueError:
            pass
    if gq and gq != '.':
        try:
            if int(gq) < 30:
                flags.append('LOW_GQ')
        except ValueError:
            pass

    flags_str = ','.join(flags) if flags else ''
    if not tier:
        continue  # Only output tiered variants

    print(f'{chrom}\t{pos}\t{ref}\t{alt}\t{gene}\t{consequence}\t{tier}\t{flags_str}\t{clinvar}\t{gnomad_af}\t{cadd}\t{spliceai}\t{loftee}\t{alphamissense}\t{revel}\t{pext}\t{zygosity}\t{dp}\t{gq}\t{ab}')
" > "${TIERED_ALL}"

# Split into tier-specific files
head -1 "${TIERED_ALL}" > "${TIER1_TSV}"
head -1 "${TIERED_ALL}" > "${TIER2_TSV}"
head -1 "${TIERED_ALL}" > "${TIER3_TSV}"

awk -F'\t' '$7 == "TIER1"' "${TIERED_ALL}" >> "${TIER1_TSV}"
awk -F'\t' '$7 == "TIER2"' "${TIERED_ALL}" >> "${TIER2_TSV}"
awk -F'\t' '$7 == "TIER3"' "${TIERED_ALL}" >> "${TIER3_TSV}"

TIER1_COUNT=$(($(wc -l < "${TIER1_TSV}" | tr -d ' ') - 1))
TIER2_COUNT=$(($(wc -l < "${TIER2_TSV}" | tr -d ' ') - 1))
TIER3_COUNT=$(($(wc -l < "${TIER3_TSV}" | tr -d ' ') - 1))
TOTAL_TIERED=$((TIER1_COUNT + TIER2_COUNT + TIER3_COUNT))

echo "  -> Tier 1 (diagnostic):  ${TIER1_COUNT} variants"
echo "  -> Tier 2 (strong):      ${TIER2_COUNT} variants"
echo "  -> Tier 3 (moderate):    ${TIER3_COUNT} variants"
echo "  -> Total tiered:         ${TOTAL_TIERED} variants"
echo

# =============================================================================
# Step 5: Kosmos artifact comparison
# =============================================================================

echo "[Step 5/8] Running Kosmos artifact check..."

if [ -x "${KOSMOS_SCRIPT}" ]; then
    bash "${KOSMOS_SCRIPT}" "${FLAGGED_VCF}" "${OUTPUT_DIR}" 2>&1 | sed 's/^/  /'
else
    echo "  [WARN] kosmos_check.sh not found or not executable. Skipping."
fi
echo

# =============================================================================
# Step 6: Oligodontia deep-dive
# =============================================================================

echo "[Step 6/8] Running oligodontia module..."

# Key oligodontia genes: WNT10A, PAX9, MSX1, EDA, EDAR, EDARADD, AXIN2, LRP6, WNT10B, LTBP3
OLIGO_GENES="WNT10A|PAX9|MSX1|EDA|EDAR|EDARADD|AXIN2|LRP6|WNT10B|LTBP3"
OLIGO_TSV="${OUTPUT_DIR}/oligodontia_variants.tsv"

head -1 "${TIERED_ALL}" > "${OLIGO_TSV}"

# Also query the full flagged VCF for any rare variants in these genes, not just tiered
${BCFTOOLS} query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CSQ_SYMBOL\t%INFO/CSQ_Consequence\t.\t.\t%INFO/CSQ_ClinVar_CLNSIG\t%INFO/gnomAD_AF\t%INFO/CSQ_CADD_PHRED\t%INFO/CSQ_SpliceAI_DS_max\t%INFO/CSQ_LoF\t%INFO/CSQ_am_class\t%INFO/CSQ_REVEL\t%INFO/CSQ_pext\t[%GT]\t[%DP]\t[%GQ]\t[%AD]\n' \
    -i "INFO/CSQ_SYMBOL ~ \"${OLIGO_GENES}\" && (INFO/gnomAD_AF_grpmax < 0.01 || INFO/gnomAD_AF_grpmax = \".\")" \
    "${FLAGGED_VCF}" 2>/dev/null \
    | awk -F'\t' 'BEGIN{OFS="\t"} {
        # Simplified: output with OLIGO tag
        gt=$17; zygosity="OTHER";
        if(gt=="1/1" || gt=="1|1") zygosity="HOM_ALT";
        else if(gt=="0/1" || gt=="0|1" || gt=="1|0") zygosity="HET";
        flags="OLIGO_GENE";
        print $1,$2,$3,$4,$5,$6,"OLIGO",flags,$9,$10,$11,$12,$13,$14,$15,$16,zygosity,$18,$19,$20
    }' >> "${OLIGO_TSV}"

OLIGO_COUNT=$(($(wc -l < "${OLIGO_TSV}" | tr -d ' ') - 1))
echo "  -> Found ${OLIGO_COUNT} rare variants in oligodontia genes"

if [ "${OLIGO_COUNT}" -gt 0 ]; then
    echo "  -> Output: ${OLIGO_TSV}"
    echo "  -> Key genes: WNT10A, PAX9, MSX1, EDA, EDAR, EDARADD, AXIN2, LRP6, WNT10B, LTBP3"
fi
echo

# =============================================================================
# Step 7: Pigmentation module
# =============================================================================

echo "[Step 7/8] Running pigmentation module..."

# Key pigmentation genes/loci
PIGMENT_GENES="MC1R|OCA2|HERC2|SLC45A2|SLC24A5|TYR|TYRP1|ASIP|IRF4|KITLG|TPCN2|BNC2"
PIGMENT_TSV="${OUTPUT_DIR}/pigmentation_variants.tsv"

head -1 "${TIERED_ALL}" > "${PIGMENT_TSV}"

${BCFTOOLS} query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CSQ_SYMBOL\t%INFO/CSQ_Consequence\t.\t.\t%INFO/CSQ_ClinVar_CLNSIG\t%INFO/gnomAD_AF\t%INFO/CSQ_CADD_PHRED\t%INFO/CSQ_SpliceAI_DS_max\t%INFO/CSQ_LoF\t%INFO/CSQ_am_class\t%INFO/CSQ_REVEL\t%INFO/CSQ_pext\t[%GT]\t[%DP]\t[%GQ]\t[%AD]\n' \
    -i "INFO/CSQ_SYMBOL ~ \"${PIGMENT_GENES}\"" \
    "${FLAGGED_VCF}" 2>/dev/null \
    | awk -F'\t' 'BEGIN{OFS="\t"} {
        gt=$17; zygosity="OTHER";
        if(gt=="1/1" || gt=="1|1") zygosity="HOM_ALT";
        else if(gt=="0/1" || gt=="0|1" || gt=="1|0") zygosity="HET";
        flags="PIGMENT_GENE";
        print $1,$2,$3,$4,$5,$6,"PIGMENT",flags,$9,$10,$11,$12,$13,$14,$15,$16,zygosity,$18,$19,$20
    }' >> "${PIGMENT_TSV}"

PIGMENT_COUNT=$(($(wc -l < "${PIGMENT_TSV}" | tr -d ' ') - 1))
echo "  -> Found ${PIGMENT_COUNT} variants in pigmentation genes"

if [ "${PIGMENT_COUNT}" -gt 0 ]; then
    echo "  -> Output: ${PIGMENT_TSV}"
fi
echo

# =============================================================================
# Summary statistics
# =============================================================================

echo "============================================================"
echo "TRIAGE SUMMARY"
echo "============================================================"

TOTAL_INPUT=$(${BCFTOOLS} view -H "${INPUT_VCF}" | wc -l | tr -d ' ')
TOTAL_FLAGGED=$(${BCFTOOLS} view -H "${FLAGGED_VCF}" | wc -l | tr -d ' ')

{
    echo "Triage Summary - $(date '+%Y-%m-%d %H:%M:%S')"
    echo "============================================================"
    echo ""
    echo "Input VCF:          ${INPUT_VCF}"
    echo "Total variants:     ${TOTAL_INPUT}"
    echo ""
    echo "--- Tiered Variants ---"
    echo "Tier 1 (diagnostic):  ${TIER1_COUNT}"
    echo "Tier 2 (strong):      ${TIER2_COUNT}"
    echo "Tier 3 (moderate):    ${TIER3_COUNT}"
    echo "Total tiered:         ${TOTAL_TIERED}"
    echo ""
    echo "--- Clustered Indels ---"
    echo "Cluster-flagged:      ${CLUSTER_COUNT}"
    echo ""
    echo "--- Module Reports ---"
    echo "Oligodontia genes:    ${OLIGO_COUNT} rare variants"
    echo "Pigmentation genes:   ${PIGMENT_COUNT} variants"
    echo ""
    echo "--- Output Files ---"
    echo "Tiered (all):    ${TIERED_ALL}"
    echo "Tier 1:          ${TIER1_TSV}"
    echo "Tier 2:          ${TIER2_TSV}"
    echo "Tier 3:          ${TIER3_TSV}"
    echo "Oligodontia:     ${OLIGO_TSV}"
    echo "Pigmentation:    ${PIGMENT_TSV}"
    echo "Kosmos report:   ${OUTPUT_DIR}/kosmos_report.txt"
    echo "Flagged VCF:     ${FLAGGED_VCF}"
} | tee "${SUMMARY}"

echo
echo "Done. Results in ${OUTPUT_DIR}/"
