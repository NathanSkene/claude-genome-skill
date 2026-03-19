#!/usr/bin/env python3
"""
HPO phenotype mapping for variant triage output.

Takes tiered variant TSVs (from local_triage.py) and enriches them with:
- HPO phenotype terms (plain English)
- Disease names (OMIM/Orphanet)
- Inheritance patterns

Usage:
    python3 hpo_annotate.py triage_output/          # annotate all tier files in directory
    python3 hpo_annotate.py tier1.tsv               # annotate single file
    python3 hpo_annotate.py --gene SCN1A            # lookup single gene
    python3 hpo_annotate.py --summary triage_output/ # grouped phenotype summary
"""

import argparse
import csv
import os
import sys
from collections import defaultdict

# Default HPO data location (relative to this script)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_HPO_FILE = os.path.join(SCRIPT_DIR, "genes_to_phenotype.txt")


def load_hpo_data(hpo_file=None):
    """Load HPO gene-to-phenotype mappings.

    Returns dict: gene_symbol -> list of {hpo_id, hpo_name, disease_id, frequency}
    """
    if hpo_file is None:
        hpo_file = DEFAULT_HPO_FILE

    if not os.path.exists(hpo_file):
        print(f"ERROR: HPO data file not found: {hpo_file}", file=sys.stderr)
        print("Download from: http://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt", file=sys.stderr)
        sys.exit(1)

    gene_phenotypes = defaultdict(list)

    with open(hpo_file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene = row["gene_symbol"]
            gene_phenotypes[gene].append({
                "hpo_id": row["hpo_id"],
                "hpo_name": row["hpo_name"],
                "disease_id": row["disease_id"],
                "frequency": row.get("frequency", "-"),
            })

    return gene_phenotypes


def extract_disease_source(disease_id):
    """Parse disease ID into source and number."""
    if ":" in disease_id:
        source, num = disease_id.split(":", 1)
        return source, num
    return "Unknown", disease_id


def get_inheritance_patterns(phenotypes):
    """Extract inheritance patterns from HPO terms."""
    inheritance_terms = {
        "HP:0000006": "AD",   # Autosomal dominant
        "HP:0000007": "AR",   # Autosomal recessive
        "HP:0001417": "XL",   # X-linked inheritance
        "HP:0001419": "XLR",  # X-linked recessive
        "HP:0001423": "XLD",  # X-linked dominant
        "HP:0001427": "MT",   # Mitochondrial
    }
    patterns = set()
    for p in phenotypes:
        if p["hpo_id"] in inheritance_terms:
            patterns.add(inheritance_terms[p["hpo_id"]])
    return patterns


def get_diseases(phenotypes):
    """Extract unique disease names from phenotype list."""
    diseases = set()
    for p in phenotypes:
        did = p["disease_id"]
        if did and did != "-":
            diseases.add(did)
    return diseases


def get_phenotype_names(phenotypes, exclude_inheritance=True):
    """Get unique phenotype names, optionally excluding inheritance terms."""
    inheritance_hpo = {
        "HP:0000006", "HP:0000007", "HP:0001417", "HP:0001419",
        "HP:0001423", "HP:0001427", "HP:0003743", "HP:0003745",
        "HP:0010982",  # Polygenic inheritance
        "HP:0012275",  # Sporadic
        "HP:0001426",  # Multifactorial inheritance
        "HP:0001428",  # Somatic mutation
        "HP:0001466",  # Contiguous gene syndrome
    }
    names = set()
    for p in phenotypes:
        if exclude_inheritance and p["hpo_id"] in inheritance_hpo:
            continue
        names.add(p["hpo_name"])
    return sorted(names)


def categorize_phenotypes(phenotype_names):
    """Group phenotypes into broad categories for summary display."""
    categories = {
        "Neurological": [],
        "Cardiac": [],
        "Metabolic": [],
        "Skeletal": [],
        "Ophthalmological": [],
        "Dermatological": [],
        "Hematological": [],
        "Renal": [],
        "Immunological": [],
        "Developmental": [],
        "Other": [],
    }

    category_keywords = {
        "Neurological": ["seizure", "epilep", "ataxia", "neuropath", "intellectual disability",
                         "spastic", "cerebr", "brain", "neurodegen", "dementia", "dystonia",
                         "tremor", "muscle weakness", "hypotonia", "migraine"],
        "Cardiac": ["cardiac", "cardiomyo", "arrhyth", "heart", "aortic", "ventricular",
                     "prolonged QT", "sudden death"],
        "Metabolic": ["metabol", "diabetes", "hyperglyc", "hypoglyc", "acidosis",
                       "aminoacid", "lipid", "cholesterol"],
        "Skeletal": ["skeletal", "scoliosis", "osteo", "fracture", "short stature",
                      "joint", "bone", "dental", "tooth", "oligodontia"],
        "Ophthalmological": ["retinal", "visual", "optic", "blindness", "cataract",
                              "glaucoma", "nystagmus", "macular"],
        "Dermatological": ["skin", "hair", "pigment", "ichthyosis", "eczema",
                            "alopecia", "nail"],
        "Hematological": ["anemia", "thrombocyt", "leukocyt", "bleeding",
                           "coagul", "hemolytic", "pancytopenia"],
        "Renal": ["renal", "kidney", "nephro", "proteinuria", "hematuria"],
        "Immunological": ["immun", "autoimmun", "infection", "lympho"],
        "Developmental": ["developmental delay", "growth", "failure to thrive",
                           "microcephaly", "macrocephaly"],
    }

    for name in phenotype_names:
        name_lower = name.lower()
        categorized = False
        for cat, keywords in category_keywords.items():
            if any(kw in name_lower for kw in keywords):
                categories[cat].append(name)
                categorized = True
                break
        if not categorized:
            categories["Other"].append(name)

    return {k: v for k, v in categories.items() if v}


def lookup_gene(gene_symbol, gene_phenotypes):
    """Pretty-print HPO data for a single gene."""
    phenotypes = gene_phenotypes.get(gene_symbol, [])
    if not phenotypes:
        print(f"No HPO phenotypes found for {gene_symbol}")
        return

    inheritance = get_inheritance_patterns(phenotypes)
    diseases = get_diseases(phenotypes)
    names = get_phenotype_names(phenotypes)
    categories = categorize_phenotypes(names)

    print(f"\n{'=' * 60}")
    print(f"Gene: {gene_symbol}")
    print(f"{'=' * 60}")
    print(f"Inheritance: {', '.join(sorted(inheritance)) if inheritance else 'Unknown'}")
    print(f"Diseases ({len(diseases)}): {', '.join(sorted(diseases)[:10])}")
    if len(diseases) > 10:
        print(f"  ... and {len(diseases) - 10} more")
    print(f"\nPhenotypes ({len(names)}):")

    for cat, phenos in categories.items():
        print(f"\n  {cat}:")
        for p in phenos[:5]:
            print(f"    - {p}")
        if len(phenos) > 5:
            print(f"    ... and {len(phenos) - 5} more")


def annotate_triage_file(input_tsv, gene_phenotypes, output_tsv=None):
    """Add HPO columns to a triage TSV file."""
    if output_tsv is None:
        base, ext = os.path.splitext(input_tsv)
        output_tsv = f"{base}.hpo{ext}"

    with open(input_tsv, "r") as fin, open(output_tsv, "w") as fout:
        reader = csv.reader(fin, delimiter="\t")
        writer = csv.writer(fout, delimiter="\t")

        header = next(reader)
        gene_col = header.index("GENE")
        writer.writerow(header + ["INHERITANCE", "DISEASES", "HPO_PHENOTYPES", "PHENOTYPE_CATEGORIES"])

        for row in reader:
            gene = row[gene_col]
            phenotypes = gene_phenotypes.get(gene, [])

            inheritance = get_inheritance_patterns(phenotypes)
            diseases = get_diseases(phenotypes)
            names = get_phenotype_names(phenotypes)
            categories = categorize_phenotypes(names)

            cat_summary = "; ".join(f"{k}({len(v)})" for k, v in categories.items())

            writer.writerow(row + [
                ",".join(sorted(inheritance)) if inheritance else "",
                ",".join(sorted(diseases)[:5]) if diseases else "",
                "; ".join(names[:10]) if names else "",
                cat_summary,
            ])

    return output_tsv


def generate_summary(triage_dir, gene_phenotypes):
    """Generate a phenotype-grouped summary across all tiers."""
    all_variants = []

    for tier_file in ["tier1.tsv", "tier2.tsv", "tier3.tsv"]:
        path = os.path.join(triage_dir, tier_file)
        if not os.path.exists(path):
            continue
        with open(path, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                row["_tier_file"] = tier_file
                all_variants.append(row)

    if not all_variants:
        print("No tiered variant files found.")
        return

    # Group by phenotype category
    category_variants = defaultdict(list)
    uncategorized = []

    for var in all_variants:
        gene = var.get("GENE", "")
        phenotypes = gene_phenotypes.get(gene, [])
        names = get_phenotype_names(phenotypes)
        categories = categorize_phenotypes(names)

        if not categories:
            uncategorized.append(var)
            continue

        for cat in categories:
            category_variants[cat].append(var)

    print("\n" + "=" * 70)
    print("PHENOTYPE-GROUPED VARIANT SUMMARY")
    print("=" * 70)

    for cat in ["Neurological", "Cardiac", "Metabolic", "Skeletal", "Ophthalmological",
                "Dermatological", "Hematological", "Renal", "Immunological",
                "Developmental", "Other"]:
        variants = category_variants.get(cat, [])
        if not variants:
            continue

        print(f"\n{'─' * 60}")
        print(f"  {cat.upper()} ({len(variants)} variants)")
        print(f"{'─' * 60}")

        for var in variants:
            gene = var.get("GENE", "?")
            csq = var.get("CONSEQUENCE", "?")
            tier = var.get("TIER", "?")
            clinvar = var.get("CLINVAR", "")
            zyg = var.get("ZYGOSITY", "?")
            hgvsp = var.get("HGVSp", "")

            phenotypes = gene_phenotypes.get(gene, [])
            names = get_phenotype_names(phenotypes)
            inheritance = get_inheritance_patterns(phenotypes)
            top_phenos = "; ".join(names[:3])

            print(f"  T{tier} {gene:15s} {csq[:30]:30s} {zyg:8s} {','.join(sorted(inheritance)):5s}")
            if clinvar:
                print(f"       ClinVar: {clinvar[:50]}")
            if hgvsp:
                print(f"       Protein: {hgvsp}")
            if top_phenos:
                print(f"       Phenotypes: {top_phenos}")
            print()

    if uncategorized:
        print(f"\n{'─' * 60}")
        print(f"  NO HPO DATA ({len(uncategorized)} variants)")
        print(f"{'─' * 60}")
        for var in uncategorized:
            gene = var.get("GENE", "?")
            csq = var.get("CONSEQUENCE", "?")
            tier = var.get("TIER", "?")
            print(f"  T{tier} {gene:15s} {csq[:30]:30s}")


def main():
    parser = argparse.ArgumentParser(description="HPO phenotype annotation for variant triage")
    parser.add_argument("input", nargs="?", help="Triage TSV file or directory")
    parser.add_argument("--gene", help="Look up phenotypes for a single gene")
    parser.add_argument("--summary", action="store_true", help="Generate phenotype-grouped summary")
    parser.add_argument("--hpo-data", default=None, help="Path to genes_to_phenotype.txt")
    args = parser.parse_args()

    gene_phenotypes = load_hpo_data(args.hpo_data)
    print(f"Loaded HPO data: {len(gene_phenotypes)} genes", file=sys.stderr)

    if args.gene:
        lookup_gene(args.gene.upper(), gene_phenotypes)
        return

    if not args.input:
        parser.print_help()
        sys.exit(1)

    if args.summary:
        if os.path.isdir(args.input):
            generate_summary(args.input, gene_phenotypes)
        else:
            print("--summary requires a directory containing tier1/2/3.tsv files")
            sys.exit(1)
        return

    # Annotate mode
    if os.path.isdir(args.input):
        for fname in ["tier1.tsv", "tier2.tsv", "tier3.tsv", "tiered_variants.tsv",
                       "oligodontia_variants.tsv", "acmg_secondary.tsv", "pigmentation_variants.tsv"]:
            path = os.path.join(args.input, fname)
            if os.path.exists(path):
                out = annotate_triage_file(path, gene_phenotypes)
                print(f"Annotated: {out}")
    else:
        out = annotate_triage_file(args.input, gene_phenotypes)
        print(f"Annotated: {out}")


if __name__ == "__main__":
    main()
