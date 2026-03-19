#!/usr/bin/env python3
"""
validation_checklist.py - Interactive post-validation tier reassignment

Reads tiered_variants.tsv and walks through each variant for manual
IGV-based validation. Records PASS/FAIL/UNCERTAIN for each variant
based on visual inspection criteria:

  1. Bidirectional strand support
  2. No flanking mismatch clusters
  3. At least 3 alt-supporting reads
  4. MAPQ > 0 for supporting reads

Outputs validated_variants.tsv with validation status column.

Usage:
    python3 validation_checklist.py tiered_variants.tsv [-o validated_variants.tsv]
    python3 validation_checklist.py tiered_variants.tsv --resume  # resume interrupted session
"""

import argparse
import csv
import json
import os
import sys
from datetime import datetime
from pathlib import Path


VALIDATION_CRITERIA = [
    ("strand_support", "Bidirectional strand support (alt reads on both + and - strands)?"),
    ("no_mismatch_clusters", "Free of flanking mismatch clusters (no artifact signature)?"),
    ("min_alt_reads", "At least 3 alt-supporting reads?"),
    ("mapq_ok", "MAPQ > 0 for alt-supporting reads?"),
]

VALID_STATUSES = {"PASS", "FAIL", "UNCERTAIN", "SKIP"}


def load_variants(tsv_path: str) -> list[dict]:
    """Load tiered variants TSV."""
    variants = []
    with open(tsv_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            variants.append(dict(row))
    return variants


def load_progress(progress_path: str) -> dict:
    """Load saved progress from interrupted session."""
    if os.path.exists(progress_path):
        with open(progress_path, "r") as f:
            return json.load(f)
    return {"validated": {}, "last_index": 0}


def save_progress(progress_path: str, progress: dict):
    """Save current progress for crash recovery."""
    with open(progress_path, "w") as f:
        json.dump(progress, f, indent=2)


def variant_key(variant: dict) -> str:
    """Generate unique key for a variant."""
    chrom = variant.get("CHROM", variant.get("chrom", variant.get("#CHROM", "")))
    pos = variant.get("POS", variant.get("pos", variant.get("position", "")))
    ref = variant.get("REF", variant.get("ref", ""))
    alt = variant.get("ALT", variant.get("alt", ""))
    return f"{chrom}:{pos}:{ref}>{alt}"


def display_variant(variant: dict, index: int, total: int):
    """Display variant information for review."""
    chrom = variant.get("CHROM", variant.get("chrom", variant.get("#CHROM", "?")))
    pos = variant.get("POS", variant.get("pos", "?"))
    ref = variant.get("REF", variant.get("ref", "?"))
    alt = variant.get("ALT", variant.get("alt", "?"))
    gene = variant.get("GENE", variant.get("gene", variant.get("Gene", "")))
    rsid = variant.get("ID", variant.get("rsid", variant.get("rsID", "")))
    tier = variant.get("TIER", variant.get("tier", variant.get("Tier", "")))
    clinvar = variant.get("CLINVAR", variant.get("clinvar", variant.get("ClinVar", "")))
    consequence = variant.get("CONSEQUENCE", variant.get("consequence", variant.get("Consequence", "")))

    print(f"\n{'='*60}")
    print(f"  Variant {index + 1}/{total}")
    print(f"{'='*60}")
    print(f"  Location:    {chrom}:{pos}")
    print(f"  Change:      {ref} > {alt}")
    if rsid and rsid != ".":
        print(f"  rsID:        {rsid}")
    if gene:
        print(f"  Gene:        {gene}")
    if tier:
        print(f"  Tier:        {tier}")
    if clinvar:
        print(f"  ClinVar:     {clinvar}")
    if consequence:
        print(f"  Consequence: {consequence}")
    print(f"{'='*60}")


def validate_variant(variant: dict) -> dict:
    """Interactive validation of a single variant."""
    results = {}

    print("\n  Validation criteria (y/n/u for yes/no/uncertain):\n")

    for criterion_id, question in VALIDATION_CRITERIA:
        while True:
            response = input(f"    {question} [y/n/u]: ").strip().lower()
            if response in ("y", "yes"):
                results[criterion_id] = True
                break
            elif response in ("n", "no"):
                results[criterion_id] = False
                break
            elif response in ("u", "uncertain", "?"):
                results[criterion_id] = None
                break
            else:
                print("    Please enter y, n, or u")

    # Determine overall status
    all_pass = all(v is True for v in results.values())
    any_fail = any(v is False for v in results.values())

    if all_pass:
        suggested = "PASS"
    elif any_fail:
        suggested = "FAIL"
    else:
        suggested = "UNCERTAIN"

    print(f"\n  Suggested status: {suggested}")

    while True:
        override = input(f"  Accept [{suggested}] or override (PASS/FAIL/UNCERTAIN/SKIP): ").strip().upper()
        if override == "":
            status = suggested
            break
        elif override in VALID_STATUSES:
            status = override
            break
        else:
            print(f"  Valid options: {', '.join(VALID_STATUSES)}")

    # Optional notes
    notes = input("  Notes (optional, press Enter to skip): ").strip()

    return {
        "status": status,
        "criteria": results,
        "notes": notes,
    }


def write_output(variants: list[dict], validations: dict, output_path: str):
    """Write validated variants TSV."""
    if not variants:
        print("No variants to write.")
        return

    # Determine columns: original + validation columns
    orig_cols = list(variants[0].keys())
    new_cols = ["VALIDATION_STATUS", "STRAND_SUPPORT", "NO_MISMATCH_CLUSTERS",
                "MIN_ALT_READS", "MAPQ_OK", "VALIDATION_NOTES"]

    all_cols = orig_cols + [c for c in new_cols if c not in orig_cols]

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=all_cols, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()

        for variant in variants:
            key = variant_key(variant)
            row = dict(variant)

            if key in validations:
                val = validations[key]
                row["VALIDATION_STATUS"] = val["status"]
                criteria = val.get("criteria", {})
                row["STRAND_SUPPORT"] = _bool_to_str(criteria.get("strand_support"))
                row["NO_MISMATCH_CLUSTERS"] = _bool_to_str(criteria.get("no_mismatch_clusters"))
                row["MIN_ALT_READS"] = _bool_to_str(criteria.get("min_alt_reads"))
                row["MAPQ_OK"] = _bool_to_str(criteria.get("mapq_ok"))
                row["VALIDATION_NOTES"] = val.get("notes", "")
            else:
                row["VALIDATION_STATUS"] = "NOT_REVIEWED"
                for col in new_cols[1:]:
                    row[col] = ""

            writer.writerow(row)

    print(f"\nWritten: {output_path}")


def _bool_to_str(val) -> str:
    if val is True:
        return "PASS"
    elif val is False:
        return "FAIL"
    elif val is None:
        return "UNCERTAIN"
    return ""


def print_summary(validations: dict, total: int):
    """Print validation summary statistics."""
    counts = {"PASS": 0, "FAIL": 0, "UNCERTAIN": 0, "SKIP": 0, "NOT_REVIEWED": 0}
    for val in validations.values():
        status = val.get("status", "NOT_REVIEWED")
        counts[status] = counts.get(status, 0) + 1

    not_reviewed = total - sum(counts.values())
    if not_reviewed > 0:
        counts["NOT_REVIEWED"] += not_reviewed

    print(f"\n{'='*60}")
    print(f"  Validation Summary")
    print(f"{'='*60}")
    print(f"  Total variants:  {total}")
    print(f"  PASS:            {counts['PASS']}")
    print(f"  FAIL:            {counts['FAIL']}")
    print(f"  UNCERTAIN:       {counts['UNCERTAIN']}")
    print(f"  SKIP:            {counts['SKIP']}")
    print(f"  NOT_REVIEWED:    {counts['NOT_REVIEWED']}")
    print(f"{'='*60}")

    if counts["FAIL"] > 0:
        print(f"\n  {counts['FAIL']} variant(s) failed validation and should be")
        print(f"  excluded or downgraded in the final report.")

    if counts["UNCERTAIN"] > 0:
        print(f"\n  {counts['UNCERTAIN']} variant(s) are uncertain and may need")
        print(f"  additional review or orthogonal confirmation.")


def main():
    parser = argparse.ArgumentParser(
        description="Interactive post-validation tier reassignment for WGS variants"
    )
    parser.add_argument("input", help="Tiered variants TSV file")
    parser.add_argument("-o", "--output", default=None,
                        help="Output validated variants TSV (default: validated_variants.tsv)")
    parser.add_argument("--resume", action="store_true",
                        help="Resume interrupted validation session")
    parser.add_argument("--tier", default=None,
                        help="Only validate variants of this tier (e.g., 1, 2)")
    parser.add_argument("--summary-only", action="store_true",
                        help="Show summary of existing validation without interactive review")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"ERROR: File not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    output_path = args.output or os.path.join(
        os.path.dirname(args.input) or ".", "validated_variants.tsv"
    )
    progress_path = args.input + ".validation_progress.json"

    # Load variants
    variants = load_variants(args.input)
    if not variants:
        print("No variants found in input file.")
        sys.exit(0)

    print(f"\nLoaded {len(variants)} variants from {args.input}")

    # Filter by tier if requested
    if args.tier:
        tier_key = None
        for k in ["TIER", "tier", "Tier"]:
            if k in variants[0]:
                tier_key = k
                break
        if tier_key:
            filtered = [v for v in variants if str(v.get(tier_key, "")) == str(args.tier)]
            print(f"Filtered to {len(filtered)} Tier {args.tier} variants")
            review_variants = filtered
        else:
            print("WARNING: No TIER column found, reviewing all variants")
            review_variants = variants
    else:
        review_variants = variants

    # Load progress
    progress = load_progress(progress_path) if args.resume else {"validated": {}, "last_index": 0}

    if args.summary_only:
        print_summary(progress["validated"], len(variants))
        sys.exit(0)

    already_done = sum(1 for v in review_variants if variant_key(v) in progress["validated"])
    remaining = len(review_variants) - already_done
    print(f"Already validated: {already_done}")
    print(f"Remaining: {remaining}")

    if remaining == 0:
        print("\nAll variants already validated.")
        write_output(variants, progress["validated"], output_path)
        print_summary(progress["validated"], len(variants))
        sys.exit(0)

    print("\nStarting validation (Ctrl+C to save and quit)\n")
    print("For each variant, review in IGV report, then answer criteria.")
    print("Commands: q=quit and save, s=skip variant\n")

    try:
        for i, variant in enumerate(review_variants):
            key = variant_key(variant)
            if key in progress["validated"]:
                continue

            display_variant(variant, i, len(review_variants))

            # Check for quit/skip
            action = input("\n  Press Enter to validate, 's' to skip, 'q' to quit: ").strip().lower()

            if action == "q":
                print("\nSaving progress...")
                save_progress(progress_path, progress)
                break
            elif action == "s":
                progress["validated"][key] = {"status": "SKIP", "criteria": {}, "notes": "Skipped"}
                save_progress(progress_path, progress)
                continue

            result = validate_variant(variant)
            progress["validated"][key] = result
            progress["last_index"] = i + 1
            save_progress(progress_path, progress)

            print(f"  -> {result['status']}")

    except KeyboardInterrupt:
        print("\n\nInterrupted. Saving progress...")
        save_progress(progress_path, progress)
        print(f"Resume with: python3 {__file__} {args.input} --resume")

    # Write output
    write_output(variants, progress["validated"], output_path)
    print_summary(progress["validated"], len(variants))

    # Clean up progress file if all done
    all_reviewed = all(
        variant_key(v) in progress["validated"] for v in review_variants
    )
    if all_reviewed and os.path.exists(progress_path):
        os.remove(progress_path)
        print(f"\nAll variants reviewed. Progress file cleaned up.")

    print(f"\nOutput: {output_path}")


if __name__ == "__main__":
    main()
