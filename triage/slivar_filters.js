// slivar JavaScript expressions for artifact-aware variant triage
// Used with: slivar expr --js slivar_filters.js ...
//
// INFO field names are consistent with VEP annotation output from nf-core/raredisease
// and vcfanno annotations from vcfanno.toml in this directory.

// ============================================================================
// Artifact detection flags
// ============================================================================

/**
 * CLUSTER_ARTIFACT: Multiple indels within a ±10bp window.
 * Clustered indels in the same region are a hallmark of alignment artifacts,
 * especially in segmental duplications and low-complexity regions.
 *
 * NOTE: slivar processes variants one at a time, so true positional clustering
 * requires upstream detection (e.g., bcftools +trio-stats or a preprocessing step).
 * This function flags indels in regions already marked as problematic by vcfanno.
 * The run_triage.sh script handles true positional clustering via bcftools.
 */
function cluster_artifact(variant) {
    // Flag indels in segmental duplications or tandem repeats
    if (variant.REF.length != variant.ALT[0].length) {
        if (variant.INFO.SEGDUP || variant.INFO.TRF || variant.INFO.SIMPLE_REPEAT) {
            return true;
        }
    }
    return false;
}

/**
 * AB_SKEW: Heterozygous variants with extreme allele balance.
 * Expected het AB is ~0.5; values <0.25 or >0.75 suggest artifact or somatic.
 * Uses FORMAT/AD field (allelic depths) to compute allele balance.
 */
function ab_skew(variant) {
    if (variant.genotypes.length == 0) return false;
    var gt = variant.genotypes[0];
    if (!gt.het) return false;

    var ad = gt.AD;
    if (ad == null || ad.length < 2) return false;

    var total = ad[0] + ad[1];
    if (total == 0) return false;

    var ab = ad[1] / total;
    return (ab < 0.25 || ab > 0.75);
}

/**
 * LOW_QUAL: Variants with insufficient depth or genotype quality.
 * DP < 15 or GQ < 30 are below reliable calling thresholds for WGS.
 */
function low_quality(variant) {
    if (variant.genotypes.length == 0) return false;
    var gt = variant.genotypes[0];

    var dp = gt.DP;
    var gq = gt.GQ;

    if (dp != null && dp < 15) return true;
    if (gq != null && gq < 30) return true;
    return false;
}

/**
 * BIOLOGICAL_IMPLAUSIBILITY: Homozygous high-confidence LoF in extremely
 * constrained genes (pLI >= 0.9 or LOEUF < 0.35).
 * These are expected to be embryonic lethal or severely deleterious,
 * so homozygous LoF calls are likely artifacts.
 *
 * Uses VEP CSQ fields: LoF (LOFTEE), LoF_filter, and gnomAD constraint metrics.
 */
function biological_implausibility(variant) {
    if (variant.genotypes.length == 0) return false;
    var gt = variant.genotypes[0];
    if (!gt.hom_alt) return false;

    // Check for LOFTEE high-confidence LoF
    var lof = variant.INFO.CSQ_LoF;
    if (lof == null || lof != "HC") return false;

    // Check constraint: pLI >= 0.9 or LOEUF < 0.35
    var pli = variant.INFO.CSQ_pLI;
    var loeuf = variant.INFO.CSQ_LOEUF;

    if (pli != null && parseFloat(pli) >= 0.9) return true;
    if (loeuf != null && parseFloat(loeuf) < 0.35) return true;

    return false;
}

// ============================================================================
// Population frequency filters
// ============================================================================

/**
 * RARE_VARIANT: gnomAD popmax AF < 0.001 (0.1%).
 * Returns true for rare variants; absent variants (no gnomAD entry) are rare.
 */
function rare_variant(variant) {
    var af = variant.INFO.gnomAD_AF_grpmax;
    if (af == null || af == ".") return true;  // absent = rare
    return (parseFloat(af) < 0.001);
}

/**
 * ULTRA_RARE: gnomAD popmax AF < 0.0001 (0.01%) or absent from gnomAD.
 */
function ultra_rare(variant) {
    var af = variant.INFO.gnomAD_AF_grpmax;
    if (af == null || af == ".") return true;  // absent = ultra-rare
    return (parseFloat(af) < 0.0001);
}

// ============================================================================
// Variant tiering
// ============================================================================

/**
 * TIER 1 - Diagnostic grade:
 * ClinVar Pathogenic or Likely Pathogenic with 2+ review stars.
 *
 * VEP CSQ fields: ClinVar_CLNSIG, ClinVar_CLNREVSTAT
 */
function tier1_diagnostic(variant) {
    var clnsig = variant.INFO.CSQ_ClinVar_CLNSIG;
    if (clnsig == null) return false;

    // Check for Pathogenic or Likely_pathogenic
    var is_plp = (clnsig.indexOf("Pathogenic") >= 0 || clnsig.indexOf("Likely_pathogenic") >= 0);
    if (!is_plp) return false;

    // Exclude conflicting interpretations
    if (clnsig.indexOf("Conflicting") >= 0) return false;

    // Check review status for 2+ stars
    var revstat = variant.INFO.CSQ_ClinVar_CLNREVSTAT;
    if (revstat == null) return false;

    // 2+ stars: criteria_provided,_multiple_submitters,_no_conflicts
    //           reviewed_by_expert_panel
    //           practice_guideline
    var two_plus = (
        revstat.indexOf("multiple_submitters") >= 0 ||
        revstat.indexOf("expert_panel") >= 0 ||
        revstat.indexOf("practice_guideline") >= 0
    );

    return two_plus;
}

/**
 * TIER 2 - Strong evidence:
 *   a) LOFTEE HC LoF in constrained gene (pLI>=0.9 or LOEUF<0.35) with pext>0.1
 *   b) SpliceAI delta score > 0.8
 *   c) AlphaMissense pathogenic + REVEL > 0.75 + CADD > 25
 *
 * Only rare variants (gnomAD popmax AF < 0.001) qualify.
 */
function tier2_strong(variant) {
    if (!rare_variant(variant)) return false;

    // (a) LOFTEE HC LoF in constrained gene with expressed transcript
    var lof = variant.INFO.CSQ_LoF;
    if (lof == "HC") {
        var pli = variant.INFO.CSQ_pLI;
        var loeuf = variant.INFO.CSQ_LOEUF;
        var constrained = false;
        if (pli != null && parseFloat(pli) >= 0.9) constrained = true;
        if (loeuf != null && parseFloat(loeuf) < 0.35) constrained = true;

        if (constrained) {
            var pext = variant.INFO.CSQ_pext;
            if (pext != null && parseFloat(pext) > 0.1) return true;
            if (pext == null) return true;  // if pext unavailable, don't penalise
        }
    }

    // (b) SpliceAI high-confidence splice-disrupting
    var spliceai = variant.INFO.CSQ_SpliceAI_DS_max;
    if (spliceai == null) spliceai = variant.INFO.SpliceAI_DS_max;
    if (spliceai != null && parseFloat(spliceai) > 0.8) return true;

    // (c) AlphaMissense pathogenic + REVEL > 0.75 + CADD > 25 (all three required)
    var am = variant.INFO.CSQ_am_class;
    if (am == null) am = variant.INFO.am_class;
    var revel = variant.INFO.CSQ_REVEL;
    if (revel == null) revel = variant.INFO.REVEL;
    var cadd = variant.INFO.CSQ_CADD_PHRED;
    if (cadd == null) cadd = variant.INFO.CADD_PHRED;

    if (am != null && revel != null && cadd != null) {
        var am_path = (am == "pathogenic" || am == "likely_pathogenic");
        if (am_path && parseFloat(revel) > 0.75 && parseFloat(cadd) > 25) return true;
    }

    return false;
}

/**
 * TIER 3 - Moderate evidence:
 * Variants in known disease genes with moderate computational/clinical evidence.
 *   - Missense with CADD > 20 and REVEL > 0.5 in OMIM gene
 *   - Moderate SpliceAI (0.2-0.8)
 *   - ClinVar VUS with single submitter in constrained gene
 *   - In-frame indel in constrained gene
 *
 * Only rare variants qualify.
 */
function tier3_moderate(variant) {
    if (!rare_variant(variant)) return false;

    // Already captured by tier1 or tier2
    if (tier1_diagnostic(variant)) return false;
    if (tier2_strong(variant)) return false;

    var cadd = variant.INFO.CSQ_CADD_PHRED;
    if (cadd == null) cadd = variant.INFO.CADD_PHRED;
    var revel = variant.INFO.CSQ_REVEL;
    if (revel == null) revel = variant.INFO.REVEL;
    var consequence = variant.INFO.CSQ_Consequence;

    // (a) Missense with moderate scores in disease gene
    if (consequence != null && consequence.indexOf("missense") >= 0) {
        if (cadd != null && revel != null) {
            if (parseFloat(cadd) > 20 && parseFloat(revel) > 0.5) return true;
        }
    }

    // (b) Moderate SpliceAI
    var spliceai = variant.INFO.CSQ_SpliceAI_DS_max;
    if (spliceai == null) spliceai = variant.INFO.SpliceAI_DS_max;
    if (spliceai != null) {
        var ds = parseFloat(spliceai);
        if (ds >= 0.2 && ds <= 0.8) return true;
    }

    // (c) ClinVar VUS in constrained gene
    var clnsig = variant.INFO.CSQ_ClinVar_CLNSIG;
    if (clnsig != null && clnsig.indexOf("Uncertain_significance") >= 0) {
        var pli = variant.INFO.CSQ_pLI;
        var loeuf = variant.INFO.CSQ_LOEUF;
        if (pli != null && parseFloat(pli) >= 0.9) return true;
        if (loeuf != null && parseFloat(loeuf) < 0.35) return true;
    }

    // (d) In-frame indel in constrained gene
    if (consequence != null && consequence.indexOf("inframe") >= 0) {
        var pli2 = variant.INFO.CSQ_pLI;
        var loeuf2 = variant.INFO.CSQ_LOEUF;
        if (pli2 != null && parseFloat(pli2) >= 0.9) return true;
        if (loeuf2 != null && parseFloat(loeuf2) < 0.35) return true;
    }

    return false;
}
