#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# Kosmos Artifact Check
# Checks for known Dante Labs / Kosmos sequencing artifacts
#
# These genes have recurrent clustered frameshift artifacts in some WGS datasets,
# caused by alignment errors in segmental duplications or low-complexity regions.
#
# Usage: ./kosmos_check.sh <vcf> <output_dir>
# =============================================================================

VCF="${1:?Usage: $0 <vcf> <output_dir>}"
OUTPUT_DIR="${2:?Usage: $0 <vcf> <output_dir>}"
BCFTOOLS="${BCFTOOLS:-bcftools}"

REPORT="${OUTPUT_DIR}/kosmos_report.txt"

# --- Known artifact loci (GRCh38 coordinates) --------------------------------
# Gene: region, description
declare -A ARTIFACT_GENES
ARTIFACT_GENES=(
    ["SON"]="chr21:33543028-33585783"
    ["SCARF2"]="chr22:20396085-20416978"
    ["AGAP3"]="chr7:151065459-151104803"
    ["PRKDC"]="chr8:47852165-48100211"
)

declare -A ARTIFACT_DESC
ARTIFACT_DESC=(
    ["SON"]="Clustered frameshifts in exon 3 (segmental duplication region)"
    ["SCARF2"]="Clustered frameshifts (low-complexity / GC-rich)"
    ["AGAP3"]="Clustered indels (segmental duplication)"
    ["PRKDC"]="Frameshift artifacts (repeat region)"
)

CLUSTER_WINDOW=50  # bp window for clustering detection

# --- Report header ------------------------------------------------------------

{
    echo "============================================================"
    echo "KOSMOS ARTIFACT CHECK"
    echo "============================================================"
    echo "VCF:     ${VCF}"
    echo "Date:    $(date '+%Y-%m-%d %H:%M:%S')"
    echo "Window:  ${CLUSTER_WINDOW}bp clustering detection"
    echo "============================================================"
    echo
} > "${REPORT}"

ARTIFACT_FOUND=0
ARTIFACT_GENES_FOUND=""

for GENE in SON SCARF2 AGAP3 PRKDC; do
    REGION="${ARTIFACT_GENES[$GENE]}"
    DESC="${ARTIFACT_DESC[$GENE]}"

    {
        echo "--- ${GENE} (${REGION}) ---"
        echo "Known artifact: ${DESC}"
        echo
    } >> "${REPORT}"

    # Extract all variants in the region
    VARIANTS_FILE="${OUTPUT_DIR}/kosmos_${GENE}_variants.tsv"

    ${BCFTOOLS} query \
        -r "${REGION}" \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%GT]\t[%DP]\t[%GQ]\t[%AD]\n' \
        "${VCF}" 2>/dev/null > "${VARIANTS_FILE}" || true

    TOTAL_VARIANTS=$(wc -l < "${VARIANTS_FILE}" | tr -d ' ')

    if [ "${TOTAL_VARIANTS}" -eq 0 ]; then
        echo "  No variants found in region." >> "${REPORT}"
        echo >> "${REPORT}"
        continue
    fi

    # Extract indels only
    INDELS_FILE="${OUTPUT_DIR}/kosmos_${GENE}_indels.tsv"
    awk -F'\t' 'length($3) != length($4)' "${VARIANTS_FILE}" > "${INDELS_FILE}" || true
    INDEL_COUNT=$(wc -l < "${INDELS_FILE}" | tr -d ' ')

    echo "  Total variants in region: ${TOTAL_VARIANTS}" >> "${REPORT}"
    echo "  Indels in region:         ${INDEL_COUNT}" >> "${REPORT}"

    if [ "${INDEL_COUNT}" -eq 0 ]; then
        echo "  Cluster assessment:       CLEAN (no indels)" >> "${REPORT}"
        echo >> "${REPORT}"
        continue
    fi

    # Check for clustering: multiple indels within CLUSTER_WINDOW bp
    CLUSTER_RESULT=$(python3 -c "
import sys

indels = []
with open('${INDELS_FILE}') as f:
    for line in f:
        parts = line.strip().split('\t')
        pos = int(parts[1])
        ref, alt = parts[2], parts[3]
        gt = parts[5] if len(parts) > 5 else '.'
        dp = parts[6] if len(parts) > 6 else '.'
        indels.append((pos, ref, alt, gt, dp))

if len(indels) < 2:
    print('NO_CLUSTER')
    sys.exit(0)

# Sort by position
indels.sort()

# Find clusters
clusters = []
current_cluster = [indels[0]]

for i in range(1, len(indels)):
    if indels[i][0] - current_cluster[-1][0] <= ${CLUSTER_WINDOW}:
        current_cluster.append(indels[i])
    else:
        if len(current_cluster) >= 2:
            clusters.append(current_cluster)
        current_cluster = [indels[i]]

if len(current_cluster) >= 2:
    clusters.append(current_cluster)

if not clusters:
    print('NO_CLUSTER')
else:
    print(f'CLUSTER_FOUND:{len(clusters)}')
    for ci, cluster in enumerate(clusters):
        span = cluster[-1][0] - cluster[0][0]
        print(f'  Cluster {ci+1}: {len(cluster)} indels spanning {span}bp ({cluster[0][0]}-{cluster[-1][0]})')
        for pos, ref, alt, gt, dp in cluster:
            change = f'{ref}>{alt}'
            if len(ref) > len(alt):
                change = f'del {ref[len(alt):]}'
            elif len(alt) > len(ref):
                change = f'ins {alt[len(ref):]}'
            print(f'    pos={pos} {change} GT={gt} DP={dp}')
" 2>/dev/null)

    if echo "${CLUSTER_RESULT}" | grep -q "CLUSTER_FOUND"; then
        ARTIFACT_FOUND=$((ARTIFACT_FOUND + 1))
        ARTIFACT_GENES_FOUND="${ARTIFACT_GENES_FOUND} ${GENE}"

        {
            echo "  Cluster assessment:       *** ARTIFACT LIKELY ***"
            echo "${CLUSTER_RESULT}" | grep -v "CLUSTER_FOUND" | grep -v "NO_CLUSTER"
            echo
            echo "  RECOMMENDATION: These ${GENE} indels are likely alignment artifacts."
            echo "  Do NOT interpret as real loss-of-function variants."
        } >> "${REPORT}"
    else
        echo "  Cluster assessment:       CLEAN (indels not clustered)" >> "${REPORT}"
    fi
    echo >> "${REPORT}"

done

# --- Summary ------------------------------------------------------------------

{
    echo "============================================================"
    echo "SUMMARY"
    echo "============================================================"
    if [ "${ARTIFACT_FOUND}" -gt 0 ]; then
        echo "ARTIFACTS DETECTED: ${ARTIFACT_FOUND} gene(s) with clustered indels"
        echo "Affected genes:${ARTIFACT_GENES_FOUND}"
        echo
        echo "ACTION: Variants in these genes should be treated as artifacts"
        echo "and excluded from clinical interpretation."
    else
        echo "NO KOSMOS ARTIFACTS DETECTED"
        echo "All checked genes (SON, SCARF2, AGAP3, PRKDC) are clean."
    fi
    echo
} >> "${REPORT}"

# Also print to stdout
cat "${REPORT}"
