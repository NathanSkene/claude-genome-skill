#!/usr/bin/env python3
"""
Genome Query Helper for Claude Skill.

Queries personal genome VCF and annotates variants
using biothings-client (MyVariant.info).

Supports both GRCh37 and GRCh38 genome builds.

Dependencies: cyvcf2, biothings-client (pip3 install cyvcf2 biothings-client)
"""

import json
import os
import subprocess
import urllib.request
import urllib.error

# --- Build configuration ---
# Set GENOME_DIR environment variable to your genome data directory
GENOME_DIR = os.getenv("GENOME_DIR", os.path.expanduser("~/genome"))

# VCF file priority: GRCh38 annotated > GRCh38 > GRCh37 rsID > GRCh37 original
VCF_GRCH38 = os.path.join(GENOME_DIR, "grch38", "genome.grch38.vcf.gz")
VCF_GRCH38_ANNOTATED = os.path.join(GENOME_DIR, "grch38", "genome.grch38.annotated.vcf.gz")
VCF_RSID = os.path.join(GENOME_DIR, "genome.rsid.vcf.gz")
VCF_ORIG = os.path.join(GENOME_DIR, "genome.vcf.gz")

# Priority: GRCh38 annotated > GRCh38 > GRCh37 rsID > GRCh37 original
if os.path.exists(VCF_GRCH38_ANNOTATED):
    VCF_PATH = VCF_GRCH38_ANNOTATED
    GENOME_BUILD = "GRCh38"
elif os.path.exists(VCF_GRCH38):
    VCF_PATH = VCF_GRCH38
    GENOME_BUILD = "GRCh38"
elif os.path.exists(VCF_RSID):
    VCF_PATH = VCF_RSID
    GENOME_BUILD = "GRCh37"
else:
    VCF_PATH = VCF_ORIG
    GENOME_BUILD = "GRCh37"

PANELS_DIR = os.path.join(os.path.dirname(__file__), "panels")
PRS_RESULTS_FILE = os.path.join(GENOME_DIR, "prs_results.json")

# Ensembl REST API — build-appropriate
ENSEMBL_REST_37 = "https://grch37.rest.ensembl.org"
ENSEMBL_REST_38 = "https://rest.ensembl.org"
ENSEMBL_REST = ENSEMBL_REST_38 if GENOME_BUILD == "GRCh38" else ENSEMBL_REST_37


def _run_bcftools(args: list[str], timeout: int = 30) -> str:
    """Run bcftools command and return stdout."""
    result = subprocess.run(
        ["bcftools"] + args,
        capture_output=True, text=True, timeout=timeout
    )
    if result.returncode != 0 and result.stderr:
        # Some bcftools commands return 1 with no output (no matches)
        if "no records" in result.stderr.lower() or not result.stdout.strip():
            return ""
    return result.stdout.strip()


def _get_biothings_client():
    """Get MyVariant.info client via biothings-client."""
    from biothings_client import get_client
    return get_client("variant")


def get_variant_by_rsid(rsid: str) -> dict | None:
    """
    Look up a variant by rsID in the VCF.

    Returns dict with chrom, pos, id, ref, alt, gt or None if not found.
    """
    if not rsid.startswith("rs"):
        rsid = f"rs{rsid}"

    output = _run_bcftools([
        "query", "-i", f'ID="{rsid}"',
        "-f", "%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%GT]\n",
        VCF_PATH
    ])

    if not output:
        return None

    # Take first match if multiple
    line = output.split("\n")[0]
    fields = line.split("\t")
    if len(fields) < 6:
        return None

    return {
        "chrom": fields[0],
        "pos": int(fields[1]),
        "id": fields[2],
        "ref": fields[3],
        "alt": fields[4],
        "gt": fields[5]
    }


def get_variant_by_pos(chrom: str, pos: int) -> dict | None:
    """
    Look up a variant by chromosome position in the VCF.

    Returns dict with chrom, pos, id, ref, alt, gt or None if not found.
    """
    chrom = _normalize_chrom(chrom)

    output = _run_bcftools([
        "view", "-H", VCF_PATH, f"{chrom}:{pos}-{pos}"
    ])

    if not output:
        return None

    line = output.split("\n")[0]
    fields = line.split("\t")
    if len(fields) < 10:
        return None

    # Parse GT from FORMAT/SAMPLE columns
    fmt = fields[8].split(":")
    sample = fields[9].split(":")
    gt_idx = fmt.index("GT") if "GT" in fmt else 0
    gt = sample[gt_idx] if gt_idx < len(sample) else "."

    return {
        "chrom": fields[0],
        "pos": int(fields[1]),
        "id": fields[2],
        "ref": fields[3],
        "alt": fields[4],
        "gt": gt
    }


def get_region(chrom: str, start: int, end: int) -> list[dict]:
    """
    Extract all variants in a genomic region.

    Returns list of variant dicts.
    """
    chrom = _normalize_chrom(chrom)

    output = _run_bcftools([
        "query",
        "-f", "%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%GT]\n",
        "-r", f"{chrom}:{start}-{end}",
        VCF_PATH
    ])

    if not output:
        return []

    variants = []
    for line in output.split("\n"):
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) >= 6:
            variants.append({
                "chrom": fields[0],
                "pos": int(fields[1]),
                "id": fields[2],
                "ref": fields[3],
                "alt": fields[4],
                "gt": fields[5]
            })

    return variants


def _normalize_chrom(chrom: str) -> str:
    """Normalize chromosome naming for the active VCF.

    GRCh37 VCF uses plain numbers (1, 2, ...).
    GRCh38 VCF may use chr-prefixed (chr1, chr2, ...) or plain numbers
    depending on reference used. Detect from VCF header.
    """
    chrom = str(chrom).replace("chr", "")
    if GENOME_BUILD == "GRCh38" and _vcf_uses_chr_prefix():
        return f"chr{chrom}"
    return chrom


_chr_prefix_cache = None

def _vcf_uses_chr_prefix() -> bool:
    """Check if VCF uses chr-prefixed chromosome names."""
    global _chr_prefix_cache
    if _chr_prefix_cache is not None:
        return _chr_prefix_cache
    try:
        output = _run_bcftools(["view", "-H", VCF_PATH, "--regions", "chr1:1-1"], timeout=5)
        _chr_prefix_cache = True
    except Exception:
        _chr_prefix_cache = False
    if not output:
        try:
            output = _run_bcftools(["view", "-H", VCF_PATH, "--regions", "1:1-1"], timeout=5)
            _chr_prefix_cache = not bool(output) if output == "" else False
        except Exception:
            _chr_prefix_cache = False
    return _chr_prefix_cache


def gene_to_coords(gene_symbol: str) -> tuple[str, int, int] | None:
    """
    Convert gene symbol to coordinates via Ensembl REST API.

    Uses GRCh37 or GRCh38 API based on active VCF build.
    Returns (chrom, start, end) or None if not found.
    """
    url = f"{ENSEMBL_REST}/lookup/symbol/homo_sapiens/{gene_symbol}?content-type=application/json"

    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode())
            chrom = _normalize_chrom(data["seq_region_name"])
            return (chrom, data["start"], data["end"])
    except (urllib.error.HTTPError, urllib.error.URLError, KeyError):
        return None


def get_gene_variants(gene_symbol: str) -> list[dict]:
    """
    Get all variants in a gene by looking up coordinates via Ensembl.
    """
    coords = gene_to_coords(gene_symbol)
    if coords is None:
        return []
    chrom, start, end = coords
    return get_region(chrom, start, end)


def annotate_variant(variant_id: str, fields: str = None) -> dict | None:
    """
    Annotate a variant using MyVariant.info via biothings-client.

    Args:
        variant_id: rsID (e.g., "rs1805007") or HGVS notation
                    (e.g., "chr16:g.89986117C>T")
        fields: Comma-separated fields to return. Default includes
                clinvar, gnomad, CADD, dbnsfp (SIFT/PolyPhen/REVEL).

    Returns annotation dict or None.

    NOTE: This makes an external API call to MyVariant.info.
    """
    if fields is None:
        fields = "clinvar,gnomad_genome,cadd,dbnsfp,dbsnp,gwascatalog"

    mv = _get_biothings_client()
    try:
        result = mv.getvariant(variant_id, fields=fields)
        return result
    except Exception:
        return None


def annotate_variants_batch(variant_ids: list[str], fields: str = None) -> list[dict]:
    """
    Batch annotate variants using MyVariant.info.

    Args:
        variant_ids: List of rsIDs or HGVS notations (max ~1000 per call).
        fields: Comma-separated fields. Default includes key annotation sources.

    Returns list of annotation dicts.

    NOTE: This makes external API calls to MyVariant.info.
    """
    if fields is None:
        fields = "clinvar,gnomad_genome,cadd,dbnsfp,dbsnp"

    mv = _get_biothings_client()
    try:
        results = mv.getvariants(variant_ids, fields=fields)
        return results if results else []
    except Exception:
        return []


def extract_panel(panel_name: str) -> list[dict]:
    """
    Extract genotypes for a predefined SNP panel.

    Uses region queries grouped by chromosome for efficiency (avoids
    per-SNP full-file scans on unannotated VCFs).

    Args:
        panel_name: Panel name (e.g., "hirisplex", "mc1r").
                    Looks for panels/{panel_name}.json.

    Returns list of dicts with panel SNP info + genotype.
    """
    panel_file = os.path.join(PANELS_DIR, f"{panel_name}.json")
    if not os.path.exists(panel_file):
        return []

    with open(panel_file) as f:
        panel = json.load(f)

    # Group SNPs by chromosome for efficient batch queries
    by_chrom = {}
    for snp in panel["snps"]:
        chrom = str(snp["chrom"])
        by_chrom.setdefault(chrom, []).append(snp)

    # Query each chromosome region once (min_pos to max_pos)
    vcf_variants = {}  # (chrom, pos) -> variant dict
    for chrom, snps in by_chrom.items():
        positions = [s["pos"] for s in snps]
        region_start = min(positions) - 1
        region_end = max(positions) + 1
        region_vars = get_region(chrom, region_start, region_end)
        for v in region_vars:
            vcf_variants[(v["chrom"], v["pos"])] = v

    # Match panel SNPs to VCF variants by position
    results = []
    for snp in panel["snps"]:
        result = {**snp}
        chrom = str(snp["chrom"])
        pos = snp["pos"]
        variant = vcf_variants.get((chrom, pos))

        if variant:
            result["gt"] = variant["gt"]
            result["found"] = True
        else:
            # No variant at this position = homozygous reference
            result["gt"] = "0/0"
            result["found"] = True  # Position exists, just hom-ref

        results.append(result)

    return results


def scan_clinvar(min_review_stars: int = 1) -> list[dict]:
    """
    Scan VCF for variants with ClinVar pathogenic/likely pathogenic annotations.

    This queries all variants with rsIDs, then batch-checks ClinVar via MyVariant.info.

    Args:
        min_review_stars: Minimum ClinVar review stars (0-4). Default 1.

    Returns list of pathogenic variant dicts with annotation.

    NOTE: Makes external API calls to MyVariant.info.
    """
    # Get all rsIDs from VCF
    output = _run_bcftools([
        "query", "-i", 'ID~"^rs"',
        "-f", "%ID\n",
        VCF_PATH
    ], timeout=120)

    if not output:
        return []

    rsids = [line.strip() for line in output.split("\n") if line.strip().startswith("rs")]

    # Batch query in chunks of 200
    mv = _get_biothings_client()
    pathogenic = []

    for i in range(0, len(rsids), 200):
        chunk = rsids[i:i+200]
        try:
            results = mv.getvariants(chunk, fields="clinvar")
        except Exception:
            continue

        for r in results:
            if not isinstance(r, dict) or "clinvar" not in r:
                continue

            clin = r["clinvar"]
            rcv = clin.get("rcv", [])
            if isinstance(rcv, dict):
                rcv = [rcv]

            for entry in rcv:
                sig = entry.get("clinical_significance", "").lower()
                if "pathogenic" in sig:
                    # Check review status
                    review = clin.get("review", {})
                    stars = review.get("review_status_stars", 0) if isinstance(review, dict) else 0

                    if stars >= min_review_stars:
                        # Get genotype
                        variant = get_variant_by_rsid(r.get("dbsnp", {}).get("rsid", r.get("_id", "")))
                        pathogenic.append({
                            "rsid": r.get("_id", ""),
                            "significance": sig,
                            "condition": entry.get("conditions", {}).get("name", "unknown"),
                            "review_stars": stars,
                            "genotype": variant["gt"] if variant else "unknown"
                        })
                        break  # One entry per variant

    return pathogenic


def scan_lof(cadd_threshold: float = 20, af_max: float = 0.01) -> list[dict]:
    """
    Scan for putative loss-of-function variants.

    Strategy:
    1. Extract variants with rsIDs from VCF
    2. Batch-query MyVariant.info for CADD, SIFT, PolyPhen, gnomAD
    3. Filter by CADD score and allele frequency
    4. Return sorted by predicted impact

    Args:
        cadd_threshold: Minimum CADD phred score (default 20 = top 1%)
        af_max: Maximum gnomAD allele frequency (default 0.01 = 1%)

    Returns list of high-impact variant dicts.

    NOTE: Makes external API calls to MyVariant.info.
    This scans ALL rsID-annotated variants, which may make thousands of API calls.
    """
    # Get non-ref genotype variants with rsIDs (skip hom-ref)
    output = _run_bcftools([
        "query", "-i", 'ID~"^rs" && GT!="0/0" && GT!="./."',
        "-f", "%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%GT]\n",
        VCF_PATH
    ], timeout=120)

    if not output:
        return []

    variants = {}
    for line in output.split("\n"):
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) >= 6 and fields[2].startswith("rs"):
            variants[fields[2]] = {
                "chrom": fields[0],
                "pos": int(fields[1]),
                "rsid": fields[2],
                "ref": fields[3],
                "alt": fields[4],
                "gt": fields[5]
            }

    rsids = list(variants.keys())
    mv = _get_biothings_client()
    high_impact = []

    # Batch query in chunks
    for i in range(0, len(rsids), 200):
        chunk = rsids[i:i+200]
        try:
            results = mv.getvariants(chunk, fields="cadd,dbnsfp,gnomad_genome,dbsnp")
        except Exception:
            continue

        for r in results:
            if not isinstance(r, dict):
                continue

            rsid = r.get("_id", "")
            if rsid not in variants:
                continue

            # Extract CADD score
            cadd = r.get("cadd", {})
            cadd_phred = cadd.get("phred") if isinstance(cadd, dict) else None

            if cadd_phred is None or cadd_phred < cadd_threshold:
                continue

            # Check gnomAD frequency
            gnomad = r.get("gnomad_genome", {})
            af = None
            if isinstance(gnomad, dict):
                af = gnomad.get("af", {})
                if isinstance(af, dict):
                    af = af.get("af")

            if af is not None and af > af_max:
                continue

            # Extract prediction scores
            dbnsfp = r.get("dbnsfp", {})
            sift_score = None
            polyphen_score = None
            revel_score = None
            gene_name = None

            if isinstance(dbnsfp, dict):
                sift = dbnsfp.get("sift", {})
                if isinstance(sift, dict):
                    sift_score = sift.get("score")
                elif isinstance(sift, list) and sift:
                    sift_score = sift[0].get("score")

                pp = dbnsfp.get("polyphen2", {})
                if isinstance(pp, dict):
                    hdiv = pp.get("hdiv", {})
                    if isinstance(hdiv, dict):
                        polyphen_score = hdiv.get("score")

                revel = dbnsfp.get("revel", {})
                if isinstance(revel, dict):
                    revel_score = revel.get("score")

                gene_name = dbnsfp.get("genename")
                if isinstance(gene_name, list):
                    gene_name = gene_name[0]

            entry = {
                **variants[rsid],
                "gene": gene_name,
                "cadd_phred": cadd_phred,
                "sift": sift_score,
                "polyphen": polyphen_score,
                "revel": revel_score,
                "gnomad_af": af
            }
            high_impact.append(entry)

    # Sort by CADD score descending
    high_impact.sort(key=lambda x: x.get("cadd_phred", 0), reverse=True)
    return high_impact


def get_gene_impact(gene_symbol: str) -> list[dict]:
    """
    Show all variants in a gene with predicted effect scores.

    1. Look up gene coordinates via Ensembl
    2. Extract all variants in the region
    3. Batch-annotate via MyVariant.info

    Returns list of annotated variant dicts sorted by CADD score.

    NOTE: Makes external API calls to Ensembl REST + MyVariant.info.
    """
    coords = gene_to_coords(gene_symbol)
    if coords is None:
        return []

    chrom, start, end = coords
    variants = get_region(chrom, start, end)

    if not variants:
        return []

    # Collect rsIDs for batch annotation
    rsids = [v["id"] for v in variants if v["id"].startswith("rs")]
    # For variants without rsIDs, use HGVS notation
    hgvs_map = {}
    for v in variants:
        if not v["id"].startswith("rs"):
            hgvs = f"chr{v['chrom']}:g.{v['pos']}{v['ref']}>{v['alt']}"
            hgvs_map[hgvs] = v

    all_ids = rsids + list(hgvs_map.keys())
    if not all_ids:
        return variants

    # Batch annotate
    mv = _get_biothings_client()
    annotated = []

    for i in range(0, len(all_ids), 200):
        chunk = all_ids[i:i+200]
        try:
            results = mv.getvariants(chunk, fields="cadd,dbnsfp,clinvar,gnomad_genome")
        except Exception:
            continue

        for r in results:
            if not isinstance(r, dict):
                continue

            rid = r.get("_id", "")

            # Find matching variant
            match = None
            for v in variants:
                if v["id"] == rid:
                    match = v
                    break
            if match is None and rid in hgvs_map:
                match = hgvs_map[rid]
            if match is None:
                continue

            # Extract annotations
            cadd = r.get("cadd", {})
            dbnsfp = r.get("dbnsfp", {})
            clinvar = r.get("clinvar", {})
            gnomad = r.get("gnomad_genome", {})

            entry = {
                **match,
                "cadd_phred": cadd.get("phred") if isinstance(cadd, dict) else None,
                "clinvar_sig": None,
                "gnomad_af": None,
                "sift": None,
                "polyphen": None,
                "gene": gene_symbol,
            }

            if isinstance(clinvar, dict):
                rcv = clinvar.get("rcv", {})
                if isinstance(rcv, list) and rcv:
                    entry["clinvar_sig"] = rcv[0].get("clinical_significance")
                elif isinstance(rcv, dict):
                    entry["clinvar_sig"] = rcv.get("clinical_significance")

            if isinstance(gnomad, dict):
                af = gnomad.get("af", {})
                if isinstance(af, dict):
                    entry["gnomad_af"] = af.get("af")
                elif isinstance(af, (int, float)):
                    entry["gnomad_af"] = af

            if isinstance(dbnsfp, dict):
                sift = dbnsfp.get("sift", {})
                if isinstance(sift, dict):
                    entry["sift"] = sift.get("score")
                pp = dbnsfp.get("polyphen2", {})
                if isinstance(pp, dict):
                    hdiv = pp.get("hdiv", {})
                    if isinstance(hdiv, dict):
                        entry["polyphen"] = hdiv.get("score")

            annotated.append(entry)

    # Sort by CADD score (None values last)
    annotated.sort(key=lambda x: x.get("cadd_phred") or 0, reverse=True)
    return annotated


def get_build_info() -> dict:
    """Return info about active genome build and VCF."""
    return {
        "build": GENOME_BUILD,
        "vcf_path": VCF_PATH,
        "ensembl_api": ENSEMBL_REST,
        "grch38_available": os.path.exists(VCF_GRCH38) or os.path.exists(VCF_GRCH38_ANNOTATED),
        "grch37_available": os.path.exists(VCF_RSID) or os.path.exists(VCF_ORIG),
    }


def get_vcf_stats() -> dict:
    """
    Get basic VCF statistics.

    Returns dict with variant count, file path, etc.
    """
    # Total variants
    total = _run_bcftools(["view", "-H", VCF_PATH, "--no-header"], timeout=120)
    total_count = len(total.split("\n")) if total else 0

    # Count with rsIDs
    rsid_output = _run_bcftools([
        "query", "-i", 'ID~"^rs"', "-f", "%ID\n", VCF_PATH
    ], timeout=120)
    rsid_count = len([l for l in rsid_output.split("\n") if l.strip()]) if rsid_output else 0

    # Stats summary
    stats_output = _run_bcftools(["stats", VCF_PATH], timeout=120)
    snp_count = None
    indel_count = None

    if stats_output:
        for line in stats_output.split("\n"):
            if line.startswith("SN\t0\tnumber of SNPs:"):
                snp_count = int(line.split("\t")[-1])
            elif line.startswith("SN\t0\tnumber of indels:"):
                indel_count = int(line.split("\t")[-1])

    return {
        "vcf_path": VCF_PATH,
        "total_variants": total_count,
        "rsid_annotated": rsid_count,
        "snps": snp_count,
        "indels": indel_count,
    }


# --- Formatted output functions for Claude ---

def format_variant(rsid: str) -> str:
    """Get formatted variant info for Claude output."""
    variant = get_variant_by_rsid(rsid)
    if variant is None:
        return f"Variant {rsid} not found in your genome VCF."

    lines = [f"## {rsid}"]
    lines.append(f"**Position:** chr{variant['chrom']}:{variant['pos']}")
    lines.append(f"**Ref/Alt:** {variant['ref']}/{variant['alt']}")
    lines.append(f"**Your genotype:** {variant['gt']}")

    gt = variant['gt']
    if gt == "0/0":
        lines.append("**Status:** Homozygous reference (you don't carry this variant)")
    elif gt == "0/1" or gt == "1/0":
        lines.append("**Status:** Heterozygous (you carry one copy)")
    elif gt == "1/1":
        lines.append("**Status:** Homozygous alternate (you carry two copies)")

    # Try annotation
    ann = annotate_variant(rsid)
    if ann:
        lines.append("")
        lines.append("### Annotation (MyVariant.info)")

        # ClinVar
        clin = ann.get("clinvar")
        if clin and isinstance(clin, dict):
            rcv = clin.get("rcv", {})
            if isinstance(rcv, list) and rcv:
                sig = rcv[0].get("clinical_significance", "unknown")
                cond = rcv[0].get("conditions", {}).get("name", "unknown")
            elif isinstance(rcv, dict):
                sig = rcv.get("clinical_significance", "unknown")
                cond = rcv.get("conditions", {}).get("name", "unknown")
            else:
                sig = "unknown"
                cond = "unknown"
            lines.append(f"- **ClinVar:** {sig} — {cond}")

        # gnomAD
        gnomad = ann.get("gnomad_genome")
        if gnomad and isinstance(gnomad, dict):
            af = gnomad.get("af", {})
            if isinstance(af, dict):
                freq = af.get("af")
            else:
                freq = af
            if freq is not None:
                lines.append(f"- **gnomAD frequency:** {freq:.4f} ({freq*100:.2f}%)")

        # CADD
        cadd = ann.get("cadd")
        if cadd and isinstance(cadd, dict):
            phred = cadd.get("phred")
            if phred is not None:
                lines.append(f"- **CADD phred:** {phred}")

        # dbSNP
        dbsnp = ann.get("dbsnp")
        if dbsnp and isinstance(dbsnp, dict):
            gene = dbsnp.get("gene", {})
            if isinstance(gene, dict):
                lines.append(f"- **Gene:** {gene.get('symbol', 'unknown')}")
            elif isinstance(gene, list) and gene:
                lines.append(f"- **Gene:** {gene[0].get('symbol', 'unknown')}")

    return "\n".join(lines)


def format_panel(panel_name: str) -> str:
    """Get formatted panel results for Claude output."""
    results = extract_panel(panel_name)
    if not results:
        return f"Panel '{panel_name}' not found."

    lines = [f"## {panel_name.upper()} Panel Results"]
    lines.append("")
    lines.append("| rsID | Gene | Pos | Ref/Alt | Your GT | Trait |")
    lines.append("|------|------|-----|---------|---------|-------|")

    for snp in results:
        gt = snp['gt']
        # Highlight non-reference genotypes
        if gt not in ("0/0", "0|0", "./."):
            gt = f"**{gt}**"
        lines.append(
            f"| {snp['rsid']} | {snp.get('gene', '')} | "
            f"chr{snp['chrom']}:{snp['pos']} | "
            f"{snp['ref']}/{snp['alt']} | {gt} | "
            f"{snp.get('trait', '')} |"
        )

    found = sum(1 for r in results if r['found'])
    lines.append(f"\n**Found:** {found}/{len(results)} SNPs in VCF")

    return "\n".join(lines)


def format_region(gene_or_coords: str) -> str:
    """Get formatted region results for Claude output."""
    # Try as gene symbol first
    if ":" not in gene_or_coords:
        gene = gene_or_coords
        coords = gene_to_coords(gene)
        if coords is None:
            return f"Gene '{gene}' not found in Ensembl."
        chrom, start, end = coords
        variants = get_region(chrom, start, end)
        header = f"## Variants in {gene} (chr{chrom}:{start}-{end})"
    else:
        # Parse coordinates
        parts = gene_or_coords.replace("chr", "").split(":")
        chrom = parts[0]
        range_parts = parts[1].split("-")
        start, end = int(range_parts[0]), int(range_parts[1])
        variants = get_region(chrom, start, end)
        header = f"## Variants in chr{chrom}:{start}-{end}"

    lines = [header]
    lines.append(f"**Total variants:** {len(variants)}")
    lines.append("")

    if variants:
        lines.append("| Pos | ID | Ref/Alt | GT |")
        lines.append("|-----|----|---------|----|")
        for v in variants:
            gt = v['gt']
            if gt not in ("0/0", "0|0", "./."):
                gt = f"**{gt}**"
            lines.append(f"| {v['pos']} | {v['id']} | {v['ref']}/{v['alt']} | {gt} |")

    return "\n".join(lines)


def format_stats() -> str:
    """Get formatted VCF stats for Claude output."""
    stats = get_vcf_stats()
    build = get_build_info()
    lines = ["## VCF Statistics"]
    lines.append(f"- **Build:** {build['build']}")
    lines.append(f"- **File:** `{stats['vcf_path']}`")
    lines.append(f"- **Total variants:** {stats['total_variants']:,}")
    lines.append(f"- **rsID annotated:** {stats['rsid_annotated']:,}")
    if stats.get('snps') is not None:
        lines.append(f"- **SNPs:** {stats['snps']:,}")
    if stats.get('indels') is not None:
        lines.append(f"- **Indels:** {stats['indels']:,}")
    if build['grch38_available'] and build['grch37_available']:
        lines.append(f"- **Note:** Both GRCh37 and GRCh38 VCFs available. Using {build['build']}.")
    return "\n".join(lines)


# --- PRS (Polygenic Risk Score) functions ---

def load_prs_results() -> dict | None:
    """Load pre-computed PRS results from JSON.

    Returns the full results dict or None if file doesn't exist.
    Local-only operation — no external API calls.
    """
    if not os.path.exists(PRS_RESULTS_FILE):
        return None
    with open(PRS_RESULTS_FILE) as f:
        return json.load(f)


def search_prs(query: str) -> list[dict]:
    """Search PRS results by trait name, EFO ID, PGS ID, or category.

    Case-insensitive substring match across trait, trait_reported,
    trait_efo, category, and pgs_id fields.

    Local-only operation — no external API calls.
    """
    data = load_prs_results()
    if data is None:
        return []

    query_lower = query.lower()
    matches = []

    for score in data["scores"]:
        searchable = " ".join([
            score.get("trait", ""),
            score.get("trait_reported", ""),
            score.get("trait_efo", ""),
            score.get("category", ""),
            score.get("pgs_id", ""),
        ]).lower()

        if query_lower in searchable:
            matches.append(score)

    # Sort by percentile descending (None values last)
    matches.sort(key=lambda x: x.get("percentile") or -999, reverse=True)
    return matches


def format_prs_summary(n: int = 20) -> str:
    """Top and bottom PRS scores ranked by percentile.

    Local-only operation — no external API calls.
    """
    data = load_prs_results()
    if data is None:
        return "No PRS results found. Run pgsc_calc on HPC first, then copy prs_results.json to $GENOME_DIR."

    scores = data["scores"]
    with_pct = [s for s in scores if s.get("percentile") is not None]
    without_pct = [s for s in scores if s.get("percentile") is None]

    if not with_pct:
        with_z = [s for s in scores if s.get("z_score") is not None]
        if with_z:
            with_z.sort(key=lambda x: x["z_score"], reverse=True)
            return _format_prs_table_z(with_z, n, data)
        return f"PRS results loaded ({len(scores)} scores) but no percentiles or z-scores available."

    with_pct.sort(key=lambda x: x["percentile"], reverse=True)

    lines = ["## Polygenic Risk Scores — Summary", ""]

    # Top N
    lines.append(f"### Top {n} (Highest Percentile)")
    lines.append("")
    lines.append("| Trait | PGS ID | Percentile | Z-score | Match Rate | Category |")
    lines.append("|-------|--------|------------|---------|------------|----------|")
    for s in with_pct[:n]:
        z = f"{s['z_score']:.2f}" if s.get("z_score") is not None else "—"
        mr = f"{s['match_rate']:.2f}" if s.get("match_rate") is not None else "—"
        lines.append(f"| {s['trait']} | {s['pgs_id']} | {s['percentile']:.1f} | {z} | {mr} | {s.get('category', '')} |")

    lines.append("")

    # Bottom N
    lines.append(f"### Bottom {n} (Lowest Percentile)")
    lines.append("")
    lines.append("| Trait | PGS ID | Percentile | Z-score | Match Rate | Category |")
    lines.append("|-------|--------|------------|---------|------------|----------|")
    for s in with_pct[-n:]:
        z = f"{s['z_score']:.2f}" if s.get("z_score") is not None else "—"
        mr = f"{s['match_rate']:.2f}" if s.get("match_rate") is not None else "—"
        lines.append(f"| {s['trait']} | {s['pgs_id']} | {s['percentile']:.1f} | {z} | {mr} | {s.get('category', '')} |")

    lines.append("")
    lines.append(f"**Total scores:** {data['n_scores']} | **Build:** {data['build']} | **Computed:** {data['date_computed']}")
    if with_pct and with_pct[0].get("ancestry_pop"):
        lines.append(f"**Ancestry reference:** {with_pct[0]['ancestry_pop']}")
    if without_pct:
        lines.append(f"**Scores without percentile:** {len(without_pct)}")
    lines.append("")
    lines.append("Use `/genome prs <trait>` to search, `/genome prs detail <PGS_ID>` for full details.")

    return "\n".join(lines)


def _format_prs_table_z(scores: list[dict], n: int, data: dict) -> str:
    """Format PRS table using z-scores when percentiles aren't available."""
    lines = ["## Polygenic Risk Scores — Summary (Z-scores)", ""]

    lines.append(f"### Top {n} (Highest Z-score)")
    lines.append("")
    lines.append("| Trait | PGS ID | Z-score | Match Rate | Category |")
    lines.append("|-------|--------|---------|------------|----------|")
    for s in scores[:n]:
        mr = f"{s['match_rate']:.2f}" if s.get("match_rate") is not None else "—"
        lines.append(f"| {s['trait']} | {s['pgs_id']} | {s['z_score']:.2f} | {mr} | {s.get('category', '')} |")

    lines.append("")
    lines.append(f"### Bottom {n} (Lowest Z-score)")
    lines.append("")
    lines.append("| Trait | PGS ID | Z-score | Match Rate | Category |")
    lines.append("|-------|--------|---------|------------|----------|")
    for s in scores[-n:]:
        mr = f"{s['match_rate']:.2f}" if s.get("match_rate") is not None else "—"
        lines.append(f"| {s['trait']} | {s['pgs_id']} | {s['z_score']:.2f} | {mr} | {s.get('category', '')} |")

    lines.append("")
    lines.append(f"**Total scores:** {data['n_scores']} | **Build:** {data['build']} | **Computed:** {data['date_computed']}")
    lines.append("*Note: Ancestry-adjusted percentiles not available. Showing raw Z-scores.*")

    return "\n".join(lines)


def format_prs_search(query: str) -> str:
    """Search PRS results and format as markdown table.

    Local-only operation — no external API calls.
    """
    matches = search_prs(query)
    if not matches:
        data = load_prs_results()
        if data is None:
            return "No PRS results found. Run pgsc_calc on HPC first."
        return f"No PRS results matching '{query}'. {data['n_scores']} scores available — try a broader search."

    lines = [f"## PRS Results — '{query}'", ""]
    lines.append(f"**Matches:** {len(matches)}")
    lines.append("")

    # Show up to 50 results
    display = matches[:50]
    has_pct = any(s.get("percentile") is not None for s in display)

    if has_pct:
        lines.append("| Trait | PGS ID | Percentile | Z-score | Match Rate | Category |")
        lines.append("|-------|--------|------------|---------|------------|----------|")
        for s in display:
            pct = f"{s['percentile']:.1f}" if s.get("percentile") is not None else "—"
            z = f"{s['z_score']:.2f}" if s.get("z_score") is not None else "—"
            mr = f"{s['match_rate']:.2f}" if s.get("match_rate") is not None else "—"
            lines.append(f"| {s['trait']} | {s['pgs_id']} | {pct} | {z} | {mr} | {s.get('category', '')} |")
    else:
        lines.append("| Trait | PGS ID | Z-score | Sum | Category |")
        lines.append("|-------|--------|---------|-----|----------|")
        for s in display:
            z = f"{s['z_score']:.2f}" if s.get("z_score") is not None else "—"
            sm = f"{s['sum']:.4f}" if s.get("sum") is not None else "—"
            lines.append(f"| {s['trait']} | {s['pgs_id']} | {z} | {sm} | {s.get('category', '')} |")

    if len(matches) > 50:
        lines.append(f"\n*Showing first 50 of {len(matches)} matches.*")

    lines.append("")
    lines.append("Use `/genome prs detail <PGS_ID>` for full details on any score.")

    return "\n".join(lines)


def format_prs_detail(pgs_id: str) -> str:
    """Detailed view of a single PRS.

    Local-only operation — no external API calls.
    """
    data = load_prs_results()
    if data is None:
        return "No PRS results found. Run pgsc_calc on HPC first."

    pgs_id_upper = pgs_id.upper()
    match = None
    for score in data["scores"]:
        if score["pgs_id"].upper() == pgs_id_upper:
            match = score
            break

    if match is None:
        return f"PGS ID '{pgs_id}' not found in results. Use `/genome prs` to see available scores."

    lines = [f"## {match['pgs_id']} — {match['trait']}", ""]

    lines.append("| Field | Value |")
    lines.append("|-------|-------|")
    lines.append(f"| **PGS ID** | {match['pgs_id']} |")
    lines.append(f"| **Trait** | {match['trait']} |")
    if match.get("trait_reported"):
        lines.append(f"| **Reported trait** | {match['trait_reported']} |")
    if match.get("trait_efo"):
        lines.append(f"| **EFO ID** | {match['trait_efo']} |")
    lines.append(f"| **Category** | {match.get('category', 'Uncategorized')} |")

    lines.append(f"| | |")
    lines.append(f"| **Raw sum** | {match.get('sum', '—')} |")
    lines.append(f"| **Variants scored** | {match.get('n_variants_scored', '—')} |")
    lines.append(f"| **Variants in model** | {match.get('n_variants_total', '—')} |")
    if match.get("match_rate") is not None:
        lines.append(f"| **Match rate** | {match['match_rate']:.1%} |")

    if match.get("percentile") is not None:
        lines.append(f"| | |")
        lines.append(f"| **Percentile** | **{match['percentile']:.1f}** |")
    if match.get("z_score") is not None:
        lines.append(f"| **Z-score** | {match['z_score']:.3f} |")
    if match.get("ancestry_pop"):
        lines.append(f"| **Ancestry ref pop** | {match['ancestry_pop']} |")

    if match.get("publication"):
        lines.append(f"| | |")
        pub_str = match["publication"]
        if match.get("pub_year"):
            pub_str += f" ({match['pub_year']})"
        lines.append(f"| **Publication** | {pub_str} |")

    lines.append("")
    lines.append(f"**Build:** {data['build']} | **Computed:** {data['date_computed']}")

    # Interpretation guidance
    lines.append("")
    if match.get("percentile") is not None:
        pct = match["percentile"]
        if pct >= 95:
            lines.append("*Interpretation: Very high — top 5% of reference population.*")
        elif pct >= 80:
            lines.append("*Interpretation: Above average — top 20% of reference population.*")
        elif pct >= 20:
            lines.append("*Interpretation: Within typical range.*")
        elif pct >= 5:
            lines.append("*Interpretation: Below average — bottom 20% of reference population.*")
        else:
            lines.append("*Interpretation: Very low — bottom 5% of reference population.*")

    return "\n".join(lines)


def main():
    """CLI interface for testing."""
    import argparse

    parser = argparse.ArgumentParser(description="Query personal genome")
    subparsers = parser.add_subparsers(dest="command")

    # variant
    p_var = subparsers.add_parser("variant", help="Look up variant by rsID")
    p_var.add_argument("rsid", help="rsID (e.g., rs1805007)")

    # region
    p_reg = subparsers.add_parser("region", help="Extract variants in region")
    p_reg.add_argument("target", help="Gene symbol or chr:start-end")

    # panel
    p_panel = subparsers.add_parser("panel", help="Extract SNP panel")
    p_panel.add_argument("name", help="Panel name (hirisplex, mc1r)")

    # stats
    subparsers.add_parser("stats", help="VCF statistics")

    # prs
    p_prs = subparsers.add_parser("prs", help="Polygenic risk scores")
    p_prs.add_argument("query", nargs="?", default=None, help="Search query or 'detail'")
    p_prs.add_argument("pgs_id", nargs="?", default=None, help="PGS ID (for detail view)")

    args = parser.parse_args()

    if args.command == "variant":
        print(format_variant(args.rsid))
    elif args.command == "region":
        print(format_region(args.target))
    elif args.command == "panel":
        print(format_panel(args.name))
    elif args.command == "stats":
        print(format_stats())
    elif args.command == "prs":
        if args.query is None:
            print(format_prs_summary())
        elif args.query.lower() == "detail" and args.pgs_id:
            print(format_prs_detail(args.pgs_id))
        else:
            print(format_prs_search(args.query))
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
