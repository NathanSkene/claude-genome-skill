#!/usr/bin/env python3
"""
Variant triage from VEP-annotated VCF.

Uses bcftools +split-vep to extract fields in one pass, then tiers variants
into Tier 1 (Diagnostic), Tier 2 (Strong candidates), Tier 3 (Moderate).
Also produces special module outputs for oligodontia, ACMG SF, and pigmentation.
"""

import argparse
import subprocess
import sys
import os
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Variant triage from VEP-annotated VCF"
    )
    parser.add_argument(
        "vcf",
        nargs="?",
        default=os.path.expanduser("~/Projects/GenomeAnalysis/results/vep/GFXC087577.vep.vcf.gz"),
        help="Path to VEP-annotated VCF (default: ~/Projects/GenomeAnalysis/results/vep/GFXC087577.vep.vcf.gz)",
    )
    parser.add_argument(
        "-o", "--outdir",
        default=None,
        help="Output directory (default: ./output relative to VCF location)",
    )
    return parser.parse_args()


# --- Config set in main() from CLI args ---
VCF = None
OUTDIR = None

# Fields to extract via split-vep (by name) + genotype fields via query
# We build a combined bcftools command that gets both VEP and sample fields
FIELDS = (
    "%CHROM\t%POS\t%REF\t%ALT\t"
    "%Consequence\t%IMPACT\t%SYMBOL\t%CLIN_SIG\t"
    "%gnomADg_NFE_AF\t%CADD_PHRED\t%SIFT\t%PolyPhen\t"
    "%HGVSp\t%CANONICAL\t%BIOTYPE\t%gnomADe_NFE_AF\t%MAX_AF"
    "[\t%GT\t%DP\t%GQ]\n"
)

# Column indices in the output
(C_CHROM, C_POS, C_REF, C_ALT,
 C_CSQ, C_IMPACT, C_SYMBOL, C_CLINSIG,
 C_GNOMAD_NFE, C_CADD, C_SIFT, C_POLYPHEN,
 C_HGVSP, C_CANONICAL, C_BIOTYPE, C_GNOMADE_NFE, C_MAXAF,
 C_GT, C_DP, C_GQ) = range(20)

# --- Gene panels ---
OLIGODONTIA_GENES = {
    "WNT10A", "PAX9", "MSX1", "EDA", "EDAR", "EDARADD",
    "AXIN2", "LRP6", "WNT10B", "LTBP3"
}

PIGMENTATION_GENES = {
    "MC1R", "OCA2", "HERC2", "SLC45A2", "SLC24A5",
    "TYR", "TYRP1", "ASIP", "IRF4"
}

# ACMG SF v3.2 genes (73 genes)
ACMG_SF_GENES = {
    "BRCA1", "BRCA2", "TP53", "STK11", "MLH1", "MSH2", "MSH6", "PMS2",
    "APC", "MUTYH", "VHL", "MEN1", "RET", "PTEN", "RB1", "SDHD", "SDHAF2",
    "SDHC", "SDHB", "TSC1", "TSC2", "WT1", "NF2", "COL3A1", "FBN1",
    "TGFBR1", "TGFBR2", "SMAD3", "ACTA2", "MYLK", "MYH11", "BMPR1A",
    "SMAD4", "KCNQ1", "KCNH2", "SCN5A", "LDLR", "APOB", "PCSK9",
    "MYH7", "MYBPC3", "TNNT2", "TNNI3", "TPM1", "ACTC1", "MYL3",
    "MYL2", "LMNA", "GLA", "FLNC", "PKP2", "DSP", "DSC2", "TMEM43",
    "DSG2", "RYR2", "CASQ2", "TRDN", "TTN", "BAG3", "DES", "RYR1",
    "CACNA1S", "ATP7B", "GAA", "OTC", "BTD", "HFE", "SERPINA1",
    "RPE65", "HBB", "HBA1", "HBA2", "PALB2", "HOXB13",
}

# --- Output header ---
HEADER = [
    "CHROM", "POS", "REF", "ALT", "GENE", "CONSEQUENCE", "IMPACT", "TIER",
    "CLINVAR", "gnomAD_NFE_AF", "CADD", "SIFT", "PolyPhen", "HGVSp",
    "ZYGOSITY", "DP", "GQ"
]

HIGH_IMPACT_CSQS = {
    "frameshift_variant", "stop_gained", "splice_donor_variant",
    "splice_acceptor_variant", "start_lost", "stop_lost",
    "transcript_ablation"
}


def parse_float(val):
    """Parse a float from VEP field, returning None for missing."""
    if not val or val == "." or val == "":
        return None
    try:
        return float(val)
    except ValueError:
        return None


def gt_to_zygosity(gt):
    """Convert GT field to human-readable zygosity."""
    if not gt or gt == ".":
        return "unknown"
    alleles = gt.replace("|", "/").split("/")
    if len(alleles) != 2:
        return gt
    if alleles[0] == alleles[1]:
        if alleles[0] == "0":
            return "hom_ref"
        return "hom_alt"
    return "het"


def is_rare(gnomad_nfe, threshold=0.01):
    """Check if variant is rare (AF below threshold or missing)."""
    if gnomad_nfe is None:
        return True
    return gnomad_nfe < threshold


def has_clinvar_pathogenic(clinsig):
    """Check if ClinVar says pathogenic/likely_pathogenic (not just conflicting)."""
    if not clinsig or clinsig == ".":
        return False
    terms = clinsig.lower().replace("_", " ").split("&")
    # Check for pathogenic/likely_pathogenic
    has_path = any(t.strip() in ("pathogenic", "likely pathogenic") for t in terms)
    # Reject if ONLY "conflicting interpretations" without explicit pathogenic
    has_benign = any("benign" in t for t in terms)
    if has_benign and not has_path:
        return False
    return has_path


def is_clinvar_benign(clinsig):
    """Check if ClinVar calls it benign/likely_benign."""
    if not clinsig or clinsig == ".":
        return False
    terms = clinsig.lower().replace("_", " ").split("&")
    return any(t.strip() in ("benign", "likely benign") for t in terms)


def tier_variant(row):
    """
    Assign tier to a variant. Returns 1, 2, 3, or None (not tiered).
    """
    csq_terms = set(row[C_CSQ].split("&")) if row[C_CSQ] != "." else set()
    impact = row[C_IMPACT]
    clinsig = row[C_CLINSIG]
    gnomad_nfe = parse_float(row[C_GNOMAD_NFE])
    cadd = parse_float(row[C_CADD])
    sift = row[C_SIFT]  # e.g. "deleterious(0.01)"
    polyphen = row[C_POLYPHEN]  # e.g. "probably_damaging(0.999)"

    # --- Tier 1: Diagnostic ---
    if has_clinvar_pathogenic(clinsig) and is_rare(gnomad_nfe, 0.05) and not is_clinvar_benign(clinsig):
        return 1

    # --- Tier 2: Strong candidates ---
    # 2a: HIGH impact + rare
    if impact == "HIGH" and is_rare(gnomad_nfe, 0.01):
        return 2

    # 2b: missense + CADD > 25 + rare
    is_missense = "missense_variant" in csq_terms
    if is_missense and cadd is not None and cadd > 25 and is_rare(gnomad_nfe, 0.01):
        return 2

    # 2c: missense + SIFT deleterious + PolyPhen probably_damaging + rare
    sift_del = sift.startswith("deleterious") if sift and sift != "." else False
    pp_dam = polyphen.startswith("probably_damaging") if polyphen and polyphen != "." else False
    if is_missense and sift_del and pp_dam and is_rare(gnomad_nfe, 0.01):
        return 2

    # --- Tier 3: Moderate ---
    if impact == "MODERATE" and is_rare(gnomad_nfe, 0.01) and cadd is not None and cadd > 15:
        return 3

    return None


def make_output_row(row, tier):
    """Build the output row list."""
    gnomad_nfe = parse_float(row[C_GNOMAD_NFE])
    return [
        row[C_CHROM],
        row[C_POS],
        row[C_REF],
        row[C_ALT],
        row[C_SYMBOL],
        row[C_CSQ],
        row[C_IMPACT],
        str(tier),
        row[C_CLINSIG] if row[C_CLINSIG] != "." else "",
        f"{gnomad_nfe:.6f}" if gnomad_nfe is not None else "",
        row[C_CADD] if row[C_CADD] != "." else "",
        row[C_SIFT] if row[C_SIFT] != "." else "",
        row[C_POLYPHEN] if row[C_POLYPHEN] != "." else "",
        row[C_HGVSP] if row[C_HGVSP] != "." else "",
        gt_to_zygosity(row[C_GT]),
        row[C_DP] if row[C_DP] != "." else "",
        row[C_GQ] if row[C_GQ] != "." else "",
    ]


def main():
    global VCF, OUTDIR
    args = parse_args()
    VCF = args.vcf
    if args.outdir:
        OUTDIR = args.outdir
    else:
        OUTDIR = os.path.join(os.path.dirname(VCF), "triage_output")
    os.makedirs(OUTDIR, exist_ok=True)

    # Build bcftools command - extract VEP fields + genotype fields in one pass
    cmd = [
        "bcftools", "+split-vep", VCF,
        "-f", FIELDS,
        "-d",       # duplicate: output one line per consequence
        "-A", "tab"  # use tab as allele separator
    ]

    print(f"Running: {' '.join(cmd)}")
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    # Accumulators
    tiers = {1: [], 2: [], 3: []}
    oligo = []
    acmg = []
    pigment = []
    seen_keys = set()  # deduplicate by chrom:pos:ref:alt:gene:tier
    total_lines = 0
    skipped_non_canonical = 0

    for line in proc.stdout:
        line = line.rstrip("\n")
        if not line:
            continue
        total_lines += 1

        fields = line.split("\t")
        if len(fields) < 20:
            continue

        # Prefer canonical transcripts to avoid duplicates
        if fields[C_CANONICAL] != "YES":
            skipped_non_canonical += 1
            continue

        # Skip non-protein-coding unless HIGH impact
        if fields[C_BIOTYPE] != "protein_coding" and fields[C_IMPACT] != "HIGH":
            continue

        gene = fields[C_SYMBOL]
        gnomad_nfe = parse_float(fields[C_GNOMAD_NFE])

        # Tier the variant
        tier = tier_variant(fields)
        if tier is not None:
            key = f"{fields[C_CHROM]}:{fields[C_POS]}:{fields[C_REF]}:{fields[C_ALT]}:{gene}:{tier}"
            if key not in seen_keys:
                seen_keys.add(key)
                out_row = make_output_row(fields, tier)
                tiers[tier].append(out_row)

        # Special panels - any rare variant in panel genes
        dedup_key_base = f"{fields[C_CHROM]}:{fields[C_POS]}:{fields[C_REF]}:{fields[C_ALT]}:{gene}"

        if gene in OLIGODONTIA_GENES and is_rare(gnomad_nfe, 0.01):
            okey = f"oligo:{dedup_key_base}"
            if okey not in seen_keys:
                seen_keys.add(okey)
                t = tier if tier else "panel_only"
                oligo.append(make_output_row(fields, t))

        if gene in PIGMENTATION_GENES and is_rare(gnomad_nfe, 0.05):
            pkey = f"pig:{dedup_key_base}"
            if pkey not in seen_keys:
                seen_keys.add(pkey)
                t = tier if tier else "panel_only"
                pigment.append(make_output_row(fields, t))

        if gene in ACMG_SF_GENES:
            clinsig = fields[C_CLINSIG]
            impact = fields[C_IMPACT]
            is_plp = has_clinvar_pathogenic(clinsig)
            is_high_rare = impact == "HIGH" and is_rare(gnomad_nfe, 0.01)
            if is_plp or is_high_rare:
                akey = f"acmg:{dedup_key_base}"
                if akey not in seen_keys:
                    seen_keys.add(akey)
                    t = tier if tier else "acmg_hit"
                    acmg.append(make_output_row(fields, t))

    proc.wait()
    stderr = proc.stderr.read()
    if proc.returncode != 0:
        print(f"bcftools error (rc={proc.returncode}):\n{stderr}", file=sys.stderr)
        sys.exit(1)

    # Write outputs
    header_line = "\t".join(HEADER)

    def write_tsv(path, rows):
        with open(path, "w") as f:
            f.write(header_line + "\n")
            for r in rows:
                f.write("\t".join(r) + "\n")

    write_tsv(os.path.join(OUTDIR, "tier1.tsv"), tiers[1])
    write_tsv(os.path.join(OUTDIR, "tier2.tsv"), tiers[2])
    write_tsv(os.path.join(OUTDIR, "tier3.tsv"), tiers[3])

    combined = []
    for t in [1, 2, 3]:
        combined.extend(tiers[t])
    write_tsv(os.path.join(OUTDIR, "tiered_variants.tsv"), combined)

    write_tsv(os.path.join(OUTDIR, "oligodontia_variants.tsv"), oligo)
    write_tsv(os.path.join(OUTDIR, "acmg_secondary.tsv"), acmg)
    write_tsv(os.path.join(OUTDIR, "pigmentation_variants.tsv"), pigment)

    # --- Summary ---
    print("\n" + "=" * 60)
    print("VARIANT TRIAGE SUMMARY")
    print("=" * 60)
    print(f"Total VEP consequence lines processed: {total_lines:,}")
    print(f"Non-canonical transcripts skipped:     {skipped_non_canonical:,}")
    print()
    print(f"  Tier 1 (Diagnostic - ClinVar P/LP):  {len(tiers[1]):,}")
    print(f"  Tier 2 (Strong - HIGH/CADD>25/SIFT): {len(tiers[2]):,}")
    print(f"  Tier 3 (Moderate - CADD>15):          {len(tiers[3]):,}")
    print(f"  Combined tiered:                      {len(combined):,}")
    print()
    print(f"  Oligodontia panel hits:               {len(oligo):,}")
    print(f"  ACMG SF v3.2 hits:                    {len(acmg):,}")
    print(f"  Pigmentation panel hits:              {len(pigment):,}")
    print()
    print(f"Output written to: {OUTDIR}/")
    print("  tier1.tsv, tier2.tsv, tier3.tsv, tiered_variants.tsv")
    print("  oligodontia_variants.tsv, acmg_secondary.tsv, pigmentation_variants.tsv")

    # Show Tier 1 hits
    if tiers[1]:
        print("\n" + "-" * 60)
        print("TIER 1 VARIANTS (ClinVar Pathogenic/Likely Pathogenic):")
        print("-" * 60)
        for r in tiers[1]:
            gene, csq, clinvar, af, zyg = r[4], r[5], r[8], r[9], r[14]
            hgvsp = r[13] if r[13] else ""
            print(f"  {gene:15s} {csq:30s} {clinvar:25s} AF={af or 'N/A':10s} {zyg:8s} {hgvsp}")

    # Show ACMG hits
    if acmg:
        print("\n" + "-" * 60)
        print("ACMG SECONDARY FINDINGS:")
        print("-" * 60)
        for r in acmg:
            gene, csq, clinvar, af, zyg = r[4], r[5], r[8], r[9], r[14]
            print(f"  {gene:15s} {csq:30s} {clinvar:25s} AF={af or 'N/A':10s} {zyg}")

    # Show oligodontia
    if oligo:
        print("\n" + "-" * 60)
        print("OLIGODONTIA PANEL:")
        print("-" * 60)
        for r in oligo:
            gene, csq, cadd, af, zyg = r[4], r[5], r[10], r[9], r[14]
            hgvsp = r[13] if r[13] else ""
            print(f"  {gene:15s} {csq:30s} CADD={cadd or 'N/A':6s} AF={af or 'N/A':10s} {zyg:8s} {hgvsp}")

    # Show pigmentation
    if pigment:
        print("\n" + "-" * 60)
        print("PIGMENTATION PANEL:")
        print("-" * 60)
        for r in pigment:
            gene, csq, cadd, af, zyg = r[4], r[5], r[10], r[9], r[14]
            hgvsp = r[13] if r[13] else ""
            print(f"  {gene:15s} {csq:30s} CADD={cadd or 'N/A':6s} AF={af or 'N/A':10s} {zyg:8s} {hgvsp}")


if __name__ == "__main__":
    main()
