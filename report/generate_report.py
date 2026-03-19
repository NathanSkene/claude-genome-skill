#!/usr/bin/env python3
"""
generate_report.py - Generate clinical HTML report from WGS re-analysis

Reads validated_variants.tsv (or tiered_variants.tsv), QC stats, and panel
results to generate a self-contained HTML report and Markdown companion.

Sections:
  1. Patient Summary + Sequencing QC
  2. Primary Findings (ClinVar P/LP)
  3. Oligodontia + Ectodermal Dysplasia Analysis
  4. Pigmentation + Red Hair Genetics
  5. High-Impact Non-ClinVar Variants
  6. Loss-of-Function Burden
  7. ACMG Secondary Findings (v3.3)
  8. Appendix: Methods, Software Versions

Usage:
    python3 generate_report.py \\
        --variants validated_variants.tsv \\
        --panels-dir ../panels/ \\
        --qc-json multiqc_data.json \\
        --sample YOUR_SAMPLE \\
        --output report.html

    # Minimal (auto-detects panels directory):
    python3 generate_report.py --variants validated_variants.tsv

    # Markdown only (no jinja2 needed):
    python3 generate_report.py --variants validated_variants.tsv --markdown-only

    # From triage directory (legacy mode):
    python3 generate_report.py --triage-dir triage/output --output report/

Dependencies: jinja2 (pip3 install jinja2) for HTML rendering.
              Markdown output works without jinja2.
"""

import argparse
import csv
import json
import os
import sys
from datetime import datetime
from pathlib import Path


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_tsv(path) -> list[dict]:
    """Load a TSV file into list of dicts."""
    path = str(path)
    if not os.path.exists(path):
        return []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader)


def load_panel(path: str) -> dict:
    """Load a gene panel JSON."""
    with open(path) as f:
        return json.load(f)


def load_qc(path: str | None) -> dict:
    """Load QC metrics from MultiQC JSON or return defaults."""
    if path and os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return {
        "mean_coverage": "N/A",
        "total_variants": "N/A",
        "alignment_rate": "N/A",
        "duplication_rate": "N/A",
    }


def get_field(variant: dict, *keys, default="") -> str:
    """Get first matching field from variant dict (case-insensitive fallback)."""
    for key in keys:
        if key in variant and variant[key]:
            return variant[key]
    lower_map = {k.lower(): k for k in variant}
    for key in keys:
        if key.lower() in lower_map:
            val = variant[lower_map[key.lower()]]
            if val:
                return val
    return default


# ---------------------------------------------------------------------------
# Variant field accessors
# ---------------------------------------------------------------------------

def get_tier(v): return get_field(v, "TIER", "tier", "Tier", default="")
def get_gene(v): return get_field(v, "GENE", "gene", "Gene", "SYMBOL", default="")
def get_chrom(v): return get_field(v, "CHROM", "chrom", "#CHROM", "Chr", default="")
def get_pos(v): return get_field(v, "POS", "pos", "Position", default="")
def get_ref(v): return get_field(v, "REF", "ref", default="")
def get_alt(v): return get_field(v, "ALT", "alt", default="")
def get_rsid(v): return get_field(v, "ID", "rsid", "rsID", default=".")
def get_genotype(v): return get_field(v, "GT", "genotype", "Genotype", default="")
def get_clinvar(v): return get_field(v, "CLINVAR", "clinvar", "ClinVar", "CLNSIG", default="")
def get_consequence(v): return get_field(v, "CONSEQUENCE", "consequence", "Consequence", "EFFECT", default="")
def get_gnomad_af(v): return get_field(v, "gnomAD_AF", "GNOMAD_AF", "AF", "gnomad_af", default="")
def get_cadd(v): return get_field(v, "CADD", "CADD_PHRED", "cadd", default="")
def get_sift(v): return get_field(v, "SIFT", "sift", default="")
def get_polyphen(v): return get_field(v, "POLYPHEN", "polyphen", "PolyPhen", default="")
def get_validation(v): return get_field(v, "VALIDATION_STATUS", "validation", default="")
def get_condition(v): return get_field(v, "CONDITION", "condition", "Disease", "CLNDN", default="")
def get_loftee(v): return get_field(v, "LOFTEE", "loftee", default="")
def get_spliceai(v): return get_field(v, "SPLICEAI", "spliceai", "SpliceAI", default="")
def get_revel(v): return get_field(v, "REVEL", "revel", default="")
def get_alphamissense(v): return get_field(v, "ALPHAMISSENSE", "alphamissense", "AlphaMissense", default="")
def get_pli(v): return get_field(v, "pLI", "PLI", "pli", default="")

def get_hgvs(v):
    hgvs = get_field(v, "HGVS", "HGVSc", "hgvs", default="")
    if hgvs:
        return hgvs
    ref, alt, pos = get_ref(v), get_alt(v), get_pos(v)
    if ref and alt:
        return f"{get_chrom(v)}:{pos} {ref}>{alt}"
    return ""


# ---------------------------------------------------------------------------
# Classification helpers
# ---------------------------------------------------------------------------

def clinsig_class(clinsig: str) -> str:
    """Map ClinVar significance to CSS class."""
    cl = clinsig.lower().replace("_", " ").replace("/", " ")
    if "pathogenic" in cl and "likely" in cl:
        return "likely-pathogenic"
    if "pathogenic" in cl:
        return "pathogenic"
    if "benign" in cl:
        return "benign"
    return "vus"


def is_clinvar_plp(v: dict) -> bool:
    """Is this variant ClinVar Pathogenic or Likely Pathogenic?"""
    cs = get_clinvar(v).lower().replace("_", " ")
    return "pathogenic" in cs and "conflicting" not in cs and "benign" not in cs


def is_lof(v: dict) -> bool:
    """Is this a loss-of-function variant?"""
    cons = get_consequence(v).lower()
    lof_terms = ["frameshift", "stop_gained", "splice_donor", "splice_acceptor",
                 "start_lost", "transcript_ablation"]
    if any(t in cons for t in lof_terms):
        return True
    # Also check LOFTEE HC
    return get_loftee(v).upper() == "HC"


def is_high_impact_non_clinvar(v: dict) -> bool:
    """High impact variant not in ClinVar."""
    cs = get_clinvar(v)
    if cs and cs != "." and "not_provided" not in cs.lower():
        return False
    cons = get_consequence(v).lower()
    high_impact = ["missense", "frameshift", "stop_gained", "splice_donor",
                   "splice_acceptor", "start_lost", "inframe"]
    if not any(t in cons for t in high_impact):
        return False
    # Require CADD >= 20 if available
    cadd = get_cadd(v)
    try:
        if float(cadd) < 20:
            return False
    except (ValueError, TypeError):
        pass
    return True


def filter_by_panel(variants: list[dict], panel: dict) -> list[dict]:
    """Filter variants to those in panel genes."""
    panel_genes = {g["symbol"].upper() for g in panel["genes"]}
    return [v for v in variants if get_gene(v).upper() in panel_genes]


def get_panel_gene_info(gene_symbol: str, panel: dict) -> dict | None:
    for g in panel["genes"]:
        if g["symbol"].upper() == gene_symbol.upper():
            return g
    return None


# ---------------------------------------------------------------------------
# Report data assembly
# ---------------------------------------------------------------------------

def assemble_report_data(
    variants: list[dict],
    panels_dir: str | None,
    qc_data: dict,
    sample_id: str,
    genome_build: str,
    pipeline_info: dict,
) -> dict:
    """Assemble all data needed for the report template."""

    data = {
        "sample_id": sample_id,
        "report_date": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "date": datetime.now().strftime("%Y-%m-%d"),
        "genome_build": genome_build,
        "total_variants": len(variants),
    }

    # QC
    data["mean_coverage"] = qc_data.get("mean_coverage", "N/A")
    data["qc"] = qc_data
    data["qc_metrics"] = [
        {"name": k, "value": str(v)} for k, v in qc_data.items()
        if k != "mean_coverage"
    ]

    # Pipeline info (defaults for Dante Labs re-analysis)
    defaults = {
        "sequencing_platform": "Dante Labs WGS (Illumina NovaSeq)",
        "alignment_method": "BWA-MEM2",
        "variant_caller": "DeepVariant",
        "annotation_tools": "VEP, vcfanno, LOFTEE, SpliceAI",
        "reference_fasta": "human_g1k_v37.fasta",
        "population_db": "gnomAD v4.0",
        "clinvar_version": "2024-12",
        "cadd_version": "1.7",
        "pipeline_name": "nf-core/raredisease",
        "pipeline_version": "2.6.0",
    }
    for k, v in defaults.items():
        data[k] = pipeline_info.get(k, v)

    # Load panels
    oligo_panel = pigment_panel = acmg_panel = None
    if panels_dir:
        for name, attr in [("oligodontia.json", "oligo"), ("pigmentation.json", "pigment"), ("acmg_sf_v3.3.json", "acmg")]:
            p = os.path.join(panels_dir, name)
            if os.path.exists(p):
                panel = load_panel(p)
                if attr == "oligo": oligo_panel = panel
                elif attr == "pigment": pigment_panel = panel
                elif attr == "acmg": acmg_panel = panel

    # Tier counts
    data["tier_counts"] = {
        "tier1": sum(1 for v in variants if get_tier(v) == "1"),
        "tier2": sum(1 for v in variants if get_tier(v) == "2"),
        "tier3": sum(1 for v in variants if get_tier(v) == "3"),
        "artifacts": 0,
    }
    data["tier1_count"] = data["tier_counts"]["tier1"]
    data["tier2_count"] = data["tier_counts"]["tier2"]

    # Section 2: Primary findings (ClinVar P/LP)
    primary = []
    for v in variants:
        if is_clinvar_plp(v):
            primary.append({
                "gene": get_gene(v), "hgvs": get_hgvs(v), "rsid": get_rsid(v),
                "clinsig": get_clinvar(v), "clinsig_class": clinsig_class(get_clinvar(v)),
                "condition": get_condition(v), "genotype": get_genotype(v),
                "tier": get_tier(v) or "?",
                "validation": get_validation(v) or "NOT_REVIEWED",
            })
    data["primary_findings"] = primary
    data["tier1"] = [v for v in variants if get_tier(v) == "1"]

    # Section 3: Oligodontia
    if oligo_panel:
        oligo_vars = filter_by_panel(variants, oligo_panel)
        data["oligodontia_variants"] = data["oligodontia"] = [{
            "gene": get_gene(v), "hgvs": get_hgvs(v), "rsid": get_rsid(v),
            "consequence": get_consequence(v), "clinsig": get_clinvar(v),
            "clinsig_class": clinsig_class(get_clinvar(v)),
            "gnomad_af": get_gnomad_af(v), "genotype": get_genotype(v),
            "inheritance": (get_panel_gene_info(get_gene(v), oligo_panel) or {}).get("inheritance", ""),
            "conditions": "; ".join((get_panel_gene_info(get_gene(v), oligo_panel) or {}).get("conditions", [])),
            # Legacy fields
            "GENE": get_gene(v), "CHROM": get_chrom(v), "POS": get_pos(v),
            "REF": get_ref(v), "ALT": get_alt(v), "CONSEQUENCE": get_consequence(v),
            "CLINVAR": get_clinvar(v), "CADD": get_cadd(v), "GT": get_genotype(v),
        } for v in oligo_vars]
        data["oligodontia_monitoring"] = (
            "<p>Variants detected in oligodontia/ectodermal dysplasia genes. "
            "Consider dental panoramic radiograph if not already performed. "
            "Discuss with clinical genetics if pathogenic/likely pathogenic variants confirmed.</p>"
        ) if data["oligodontia_variants"] else ""
    else:
        data["oligodontia_variants"] = data["oligodontia"] = []
        data["oligodontia_monitoring"] = ""

    # Section 4: Pigmentation
    if pigment_panel:
        pigment_vars = filter_by_panel(variants, pigment_panel)
        data["pigmentation_variants"] = data["pigmentation"] = [{
            "gene": get_gene(v), "rsid": get_rsid(v), "hgvs": get_hgvs(v),
            "genotype": get_genotype(v),
            "effect": (get_panel_gene_info(get_gene(v), pigment_panel) or {}).get("notes", ""),
            # Legacy
            "GENE": get_gene(v), "CHROM": get_chrom(v), "POS": get_pos(v),
            "REF": get_ref(v), "ALT": get_alt(v), "GT": get_genotype(v),
            "CADD": get_cadd(v), "CONSEQUENCE": get_consequence(v),
        } for v in pigment_vars]
        data["pigmentation_summary"] = ""
    else:
        data["pigmentation_variants"] = data["pigmentation"] = []
        data["pigmentation_summary"] = ""

    # Section 5: High-impact non-ClinVar
    data["highimpact_variants"] = data["tier2"] = [{
        "gene": get_gene(v), "hgvs": get_hgvs(v),
        "consequence": get_consequence(v), "cadd": get_cadd(v),
        "sift": get_sift(v), "polyphen": get_polyphen(v),
        "gnomad_af": get_gnomad_af(v), "genotype": get_genotype(v),
        "tier": get_tier(v) or "?",
        # Legacy
        "GENE": get_gene(v), "CHROM": get_chrom(v), "POS": get_pos(v),
        "CONSEQUENCE": get_consequence(v), "LOFTEE": get_loftee(v),
        "SPLICEAI": get_spliceai(v), "ALPHAMISSENSE": get_alphamissense(v),
        "CADD": get_cadd(v), "REVEL": get_revel(v),
        "GNOMAD_AF": get_gnomad_af(v), "GT": get_genotype(v),
    } for v in variants if is_high_impact_non_clinvar(v)]

    # Section 6: Loss-of-function
    lof_list = [{
        "gene": get_gene(v), "hgvs": get_hgvs(v),
        "lof_type": get_consequence(v), "gnomad_af": get_gnomad_af(v),
        "pli": get_pli(v), "genotype": get_genotype(v),
        "validation": get_validation(v) or "NOT_REVIEWED",
        # Legacy
        "GENE": get_gene(v), "CONSEQUENCE": get_consequence(v),
        "CADD": get_cadd(v), "GNOMAD_AF": get_gnomad_af(v),
        "GT": get_genotype(v), "LOFTEE": get_loftee(v),
    } for v in variants if is_lof(v)]
    data["lof_variants"] = lof_list
    data["lof_total"] = len(lof_list)
    try:
        data["lof_constrained"] = sum(1 for v in lof_list if float(v.get("pli") or 0) > 0.9)
    except (ValueError, TypeError):
        data["lof_constrained"] = 0

    # Section 7: ACMG secondary findings
    if acmg_panel:
        acmg_genes = {g["symbol"].upper() for g in acmg_panel["genes"]}
        acmg_vars = []
        for v in variants:
            if get_gene(v).upper() in acmg_genes and is_clinvar_plp(v):
                gi = get_panel_gene_info(get_gene(v), acmg_panel) or {}
                acmg_vars.append({
                    "gene": get_gene(v), "hgvs": get_hgvs(v),
                    "clinsig": get_clinvar(v), "clinsig_class": clinsig_class(get_clinvar(v)),
                    "condition": get_condition(v) or "; ".join(gi.get("conditions", [])),
                    "inheritance": gi.get("inheritance", ""),
                    "genotype": get_genotype(v),
                    "validation": get_validation(v) or "NOT_REVIEWED",
                    # Legacy
                    "GENE": get_gene(v), "CHROM": get_chrom(v), "POS": get_pos(v),
                    "CLINVAR": get_clinvar(v), "CONSEQUENCE": get_consequence(v),
                    "GT": get_genotype(v),
                })
        data["acmg_findings"] = acmg_vars
        data["acmg_monitoring"] = (
            "<p>ACMG secondary findings detected. These are medically actionable and "
            "should be confirmed by a CLIA/UKAS-accredited clinical laboratory before "
            "any clinical action is taken.</p>"
        ) if acmg_vars else ""
    else:
        data["acmg_findings"] = []
        data["acmg_monitoring"] = ""

    return data


# ---------------------------------------------------------------------------
# HTML rendering
# ---------------------------------------------------------------------------

def render_template(template_path: str, data: dict) -> str:
    """Render HTML template with Jinja2 if available, else basic fallback."""
    try:
        from jinja2 import Environment, FileSystemLoader
        env = Environment(
            loader=FileSystemLoader(os.path.dirname(template_path)),
            autoescape=True,
        )
        template = env.get_template(os.path.basename(template_path))
        return template.render(**data)
    except ImportError:
        print("WARNING: jinja2 not installed. Using basic HTML fallback.", file=sys.stderr)
        print("Install with: pip3 install jinja2", file=sys.stderr)
        return generate_basic_html(data)


def generate_basic_html(data: dict) -> str:
    """Generate basic HTML without jinja2 (self-contained fallback)."""
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>WGS Report - {data['sample_id']}</title>
<style>
    body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; max-width: 1200px; margin: 0 auto; padding: 20px; color: #333; }}
    h1 {{ color: #1a1a2e; border-bottom: 3px solid #333; padding-bottom: 10px; }}
    h2 {{ color: #1a1a2e; border-bottom: 1px solid #ddd; padding-bottom: 5px; margin-top: 40px; }}
    .disclaimer {{ background: #fff3cd; border: 2px solid #ffc107; padding: 16px; border-radius: 8px; margin: 20px 0; }}
    .disclaimer strong {{ color: #856404; }}
    table {{ border-collapse: collapse; width: 100%; margin: 10px 0; font-size: 0.85rem; }}
    th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
    th {{ background: #f8f9fa; }}
    tr:nth-child(even) {{ background: #f9f9f9; }}
    .summary-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 12px; margin: 20px 0; }}
    .summary-card {{ background: #f8f9fa; border: 1px solid #dee2e6; border-radius: 8px; padding: 16px; text-align: center; }}
    .summary-card .value {{ font-size: 1.4rem; font-weight: 600; }}
    .summary-card .label {{ font-size: 0.8rem; color: #6c757d; text-transform: uppercase; }}
    .monitoring {{ background: #e8f4fd; border-left: 4px solid #0d6efd; padding: 12px 16px; margin: 12px 0; border-radius: 0 4px 4px 0; }}
    @media print {{ body {{ max-width: none; padding: 0; font-size: 10pt; }} table {{ font-size: 8pt; }} }}
</style>
</head>
<body>

<div class="disclaimer">
    <strong>RESEARCH USE ONLY</strong> &mdash; This report is generated from consumer-grade whole genome sequencing
    using a research re-analysis pipeline. Not clinically validated. Do not make medical decisions without clinical-grade confirmation.
</div>

<h1>Whole Genome Sequencing Analysis Report</h1>
<p><strong>Sample:</strong> {data['sample_id']} | <strong>Build:</strong> {data['genome_build']} | <strong>Date:</strong> {data['report_date']}</p>

<div class="summary-grid">
    <div class="summary-card"><div class="value">{data['mean_coverage']}</div><div class="label">Mean Coverage</div></div>
    <div class="summary-card"><div class="value">{data['total_variants']}</div><div class="label">Total Variants</div></div>
    <div class="summary-card"><div class="value">{data['tier1_count']}</div><div class="label">Tier 1</div></div>
    <div class="summary-card"><div class="value">{data['tier2_count']}</div><div class="label">Tier 2</div></div>
</div>
"""

    def make_table(rows, columns):
        if not rows:
            return "<p>None detected.</p>"
        t = "<table><thead><tr>"
        for col_name, _ in columns:
            t += f"<th>{col_name}</th>"
        t += "</tr></thead><tbody>"
        for row in rows:
            t += "<tr>"
            for _, key in columns:
                t += f"<td>{row.get(key, '')}</td>"
            t += "</tr>"
        t += "</tbody></table>"
        return t

    # Primary Findings
    html += "<h2>Primary Findings (ClinVar P/LP)</h2>"
    html += make_table(data["primary_findings"], [
        ("Gene", "gene"), ("Variant", "hgvs"), ("rsID", "rsid"),
        ("ClinVar", "clinsig"), ("Condition", "condition"),
        ("Genotype", "genotype"), ("Tier", "tier"), ("Validation", "validation"),
    ])

    # Oligodontia
    html += "<h2>Oligodontia &amp; Ectodermal Dysplasia</h2>"
    html += make_table(data["oligodontia_variants"], [
        ("Gene", "gene"), ("Variant", "hgvs"), ("rsID", "rsid"),
        ("Consequence", "consequence"), ("ClinVar", "clinsig"),
        ("gnomAD AF", "gnomad_af"), ("Genotype", "genotype"),
        ("Inheritance", "inheritance"),
    ])
    if data.get("oligodontia_monitoring"):
        html += f'<div class="monitoring">{data["oligodontia_monitoring"]}</div>'

    # Pigmentation
    html += "<h2>Pigmentation &amp; Red Hair Genetics</h2>"
    html += make_table(data["pigmentation_variants"], [
        ("Gene", "gene"), ("rsID", "rsid"), ("Variant", "hgvs"),
        ("Genotype", "genotype"), ("Effect", "effect"),
    ])

    # High-impact non-ClinVar
    html += "<h2>High-Impact Non-ClinVar Variants</h2>"
    html += '<p style="color:#856404;font-style:italic;">Research-grade predictions only. Require clinical validation.</p>'
    html += make_table(data["highimpact_variants"], [
        ("Gene", "gene"), ("Variant", "hgvs"), ("Consequence", "consequence"),
        ("CADD", "cadd"), ("SIFT", "sift"), ("PolyPhen", "polyphen"),
        ("gnomAD AF", "gnomad_af"), ("Genotype", "genotype"), ("Tier", "tier"),
    ])

    # LoF
    html += "<h2>Loss-of-Function Burden</h2>"
    html += f"<p>Total LoF variants: {data['lof_total']} | In constrained genes (pLI &gt; 0.9): {data['lof_constrained']}</p>"
    html += make_table(data["lof_variants"], [
        ("Gene", "gene"), ("Variant", "hgvs"), ("Type", "lof_type"),
        ("gnomAD AF", "gnomad_af"), ("pLI", "pli"),
        ("Genotype", "genotype"), ("Validation", "validation"),
    ])

    # ACMG
    html += "<h2>ACMG Secondary Findings (v3.3)</h2>"
    html += make_table(data["acmg_findings"], [
        ("Gene", "gene"), ("Variant", "hgvs"), ("ClinVar", "clinsig"),
        ("Condition", "condition"), ("Inheritance", "inheritance"),
        ("Genotype", "genotype"), ("Validation", "validation"),
    ])
    if data.get("acmg_monitoring"):
        html += f'<div class="monitoring">{data["acmg_monitoring"]}</div>'

    # Appendix
    html += f"""
<h2>Appendix: Methods &amp; Software</h2>
<table>
<tr><th>Component</th><th>Details</th></tr>
<tr><td>Sequencing</td><td>{data['sequencing_platform']}</td></tr>
<tr><td>Alignment</td><td>{data['alignment_method']}</td></tr>
<tr><td>Variant Calling</td><td>{data['variant_caller']}</td></tr>
<tr><td>Annotation</td><td>{data['annotation_tools']}</td></tr>
<tr><td>Reference</td><td>{data['genome_build']} ({data['reference_fasta']})</td></tr>
<tr><td>Population DB</td><td>{data['population_db']}</td></tr>
<tr><td>ClinVar</td><td>{data['clinvar_version']}</td></tr>
<tr><td>CADD</td><td>v{data['cadd_version']}</td></tr>
<tr><td>Pipeline</td><td>{data['pipeline_name']} v{data['pipeline_version']}</td></tr>
</table>
<p><em>Report generated: {data['report_date']}</em></p>
</body></html>"""

    return html


# ---------------------------------------------------------------------------
# Markdown report
# ---------------------------------------------------------------------------

def generate_markdown(data: dict) -> str:
    """Generate Markdown version of the report."""
    lines = [
        f"# WGS Analysis Report - {data['sample_id']}",
        "",
        f"**Build:** {data['genome_build']} | **Date:** {data['report_date']}",
        "",
        "> **RESEARCH USE ONLY** - This report is generated from consumer-grade WGS",
        "> using a research re-analysis pipeline. Not clinically validated.",
        "",
        "## Summary",
        "",
        f"| Metric | Value |",
        f"|--------|-------|",
        f"| Mean Coverage | {data['mean_coverage']} |",
        f"| Total Variants | {data['total_variants']} |",
        f"| Tier 1 Findings | {data['tier1_count']} |",
        f"| Tier 2 Findings | {data['tier2_count']} |",
        "",
    ]

    def add_table(title, rows, columns):
        lines.append(f"## {title}")
        lines.append("")
        if not rows:
            lines.append("None detected.")
            lines.append("")
            return
        header = "| " + " | ".join(c[0] for c in columns) + " |"
        sep = "|" + "|".join("------" for _ in columns) + "|"
        lines.append(header)
        lines.append(sep)
        for r in rows:
            row = "| " + " | ".join(str(r.get(c[1], "")) for c in columns) + " |"
            lines.append(row)
        lines.append("")

    add_table("Primary Findings (ClinVar P/LP)", data["primary_findings"], [
        ("Gene", "gene"), ("Variant", "hgvs"), ("rsID", "rsid"),
        ("ClinVar", "clinsig"), ("Condition", "condition"),
        ("Genotype", "genotype"), ("Tier", "tier"), ("Validation", "validation"),
    ])

    add_table("Oligodontia & Ectodermal Dysplasia", data["oligodontia_variants"], [
        ("Gene", "gene"), ("Variant", "hgvs"), ("Consequence", "consequence"),
        ("ClinVar", "clinsig"), ("gnomAD AF", "gnomad_af"),
        ("Genotype", "genotype"), ("Inheritance", "inheritance"),
    ])

    add_table("Pigmentation & Red Hair Genetics", data["pigmentation_variants"], [
        ("Gene", "gene"), ("rsID", "rsid"), ("Variant", "hgvs"),
        ("Genotype", "genotype"), ("Effect", "effect"),
    ])

    add_table("High-Impact Non-ClinVar Variants", data["highimpact_variants"], [
        ("Gene", "gene"), ("Variant", "hgvs"), ("Consequence", "consequence"),
        ("CADD", "cadd"), ("gnomAD AF", "gnomad_af"),
        ("Genotype", "genotype"), ("Tier", "tier"),
    ])

    lines.append("## Loss-of-Function Burden")
    lines.append("")
    lines.append(f"Total LoF: {data['lof_total']} | Constrained (pLI > 0.9): {data['lof_constrained']}")
    lines.append("")
    add_table("LoF Variants", data["lof_variants"], [
        ("Gene", "gene"), ("Variant", "hgvs"), ("Type", "lof_type"),
        ("gnomAD AF", "gnomad_af"), ("pLI", "pli"),
        ("Genotype", "genotype"), ("Validation", "validation"),
    ])

    add_table("ACMG Secondary Findings (v3.3)", data["acmg_findings"], [
        ("Gene", "gene"), ("Variant", "hgvs"), ("ClinVar", "clinsig"),
        ("Condition", "condition"), ("Inheritance", "inheritance"),
        ("Genotype", "genotype"), ("Validation", "validation"),
    ])

    lines.extend([
        "## Methods",
        "",
        f"- **Pipeline:** {data['pipeline_name']} v{data['pipeline_version']}",
        f"- **Alignment:** {data['alignment_method']}",
        f"- **Variant Calling:** {data['variant_caller']}",
        f"- **Annotation:** {data['annotation_tools']}",
        f"- **Reference:** {data['genome_build']}",
        f"- **Population DB:** {data['population_db']}",
        f"- **ClinVar:** {data['clinvar_version']}",
        "",
        f"---",
        f"*Report generated: {data['report_date']}*",
    ])

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Generate WGS clinical report (HTML + Markdown)")
    # New interface
    parser.add_argument("--variants", default=None, help="Validated/tiered variants TSV")
    parser.add_argument("--panels-dir", default=None, help="Panels directory")
    parser.add_argument("--qc-json", default=None, help="MultiQC data JSON")
    parser.add_argument("--sample", default="YOUR_SAMPLE", help="Sample ID")
    parser.add_argument("--build", default="GRCh37/hg19", help="Genome build")
    parser.add_argument("--output", default="report.html", help="Output HTML file")
    parser.add_argument("--template", default=None, help="HTML template (Jinja2)")
    parser.add_argument("--pipeline-info", default=None, help="Pipeline info JSON")
    parser.add_argument("--markdown-only", action="store_true", help="Markdown output only")
    # Legacy interface
    parser.add_argument("--triage-dir", default=None, help="Triage output directory (legacy)")
    parser.add_argument("--multiqc", default=None, help="MultiQC JSON (legacy, use --qc-json)")
    parser.add_argument("--sample-id", default=None, help="Sample ID (legacy, use --sample)")
    args = parser.parse_args()

    # Resolve legacy arguments
    sample_id = args.sample_id or args.sample
    qc_json = args.qc_json or args.multiqc

    project_dir = Path(__file__).parent.parent

    # Resolve panels directory
    panels_dir = args.panels_dir
    if panels_dir is None:
        default_panels = project_dir / "panels"
        if default_panels.exists():
            panels_dir = str(default_panels)

    # Load variants
    if args.variants:
        variants = load_tsv(args.variants)
    elif args.triage_dir:
        # Legacy: load from triage directory tier files
        triage = Path(args.triage_dir)
        variants = []
        for tier_file in ["tier1.tsv", "tier2.tsv", "tier3.tsv"]:
            tier_variants = load_tsv(triage / tier_file)
            tier_num = tier_file.replace("tier", "").replace(".tsv", "")
            for v in tier_variants:
                v["TIER"] = tier_num
            variants.extend(tier_variants)
    else:
        print("ERROR: provide --variants or --triage-dir", file=sys.stderr)
        sys.exit(1)

    print(f"Loaded {len(variants)} variants")
    if panels_dir:
        print(f"Panels: {panels_dir}")

    # Load QC and pipeline info
    qc_data = load_qc(qc_json)
    pipeline_info = {}
    if args.pipeline_info and os.path.exists(args.pipeline_info):
        with open(args.pipeline_info) as f:
            pipeline_info = json.load(f)

    # Assemble data
    data = assemble_report_data(
        variants=variants,
        panels_dir=panels_dir,
        qc_data=qc_data,
        sample_id=sample_id,
        genome_build=args.build,
        pipeline_info=pipeline_info,
    )

    # Generate Markdown
    output_path = args.output
    if output_path.endswith(".html"):
        md_path = output_path.replace(".html", ".md")
    else:
        md_path = output_path + ".md"

    md_content = generate_markdown(data)
    os.makedirs(os.path.dirname(md_path) or ".", exist_ok=True)
    with open(md_path, "w") as f:
        f.write(md_content)
    print(f"Markdown report: {md_path}")

    if args.markdown_only:
        return

    # Generate HTML
    template_path = args.template
    if template_path is None:
        template_path = os.path.join(os.path.dirname(__file__), "report_template.html")

    if os.path.exists(template_path):
        html_content = render_template(template_path, data)
    else:
        print(f"Template not found at {template_path}, using basic HTML.", file=sys.stderr)
        html_content = generate_basic_html(data)

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w") as f:
        f.write(html_content)

    size_kb = os.path.getsize(output_path) / 1024
    print(f"HTML report: {output_path} ({size_kb:.0f} KB)")
    print(f"\nReport sections:")
    print(f"  Primary findings (ClinVar P/LP): {len(data['primary_findings'])}")
    print(f"  Oligodontia panel:               {len(data['oligodontia_variants'])}")
    print(f"  Pigmentation panel:              {len(data['pigmentation_variants'])}")
    print(f"  High-impact non-ClinVar:         {len(data['highimpact_variants'])}")
    print(f"  Loss-of-function:                {data['lof_total']}")
    print(f"  ACMG secondary findings:         {len(data['acmg_findings'])}")


if __name__ == "__main__":
    main()
