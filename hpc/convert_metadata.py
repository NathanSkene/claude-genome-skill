#!/usr/bin/env python3
"""Convert bulk PGS Catalog CSV metadata to JSON format for process_prs_results.py."""

import csv
import json
import sys

scores_csv = sys.argv[1] if len(sys.argv) > 1 else "pgs_all_metadata_scores.csv"
pubs_csv = sys.argv[2] if len(sys.argv) > 2 else "pgs_all_metadata_publications.csv"
output = sys.argv[3] if len(sys.argv) > 3 else "pgs_metadata.json"

# Load publications
pubs = {}
with open(pubs_csv) as f:
    reader = csv.DictReader(f)
    for row in reader:
        pgp_id = row.get("PGS Publication/Study (PGP) ID", "")
        pub_date = row.get("Publication Date", "") or ""
        pubs[pgp_id] = {
            "author": row.get("First Author", ""),
            "year": pub_date[:4] if pub_date else "",
        }

# Load scores metadata
metadata = {}
with open(scores_csv) as f:
    reader = csv.DictReader(f)
    for row in reader:
        pgs_id = row["Polygenic Score (PGS) ID"]
        efo_labels = row.get("Mapped Trait(s) (EFO label)", "") or ""
        efo_ids = row.get("Mapped Trait(s) (EFO ID)", "") or ""
        reported = row.get("Reported Trait", "") or ""
        n_var_str = row.get("Number of Variants", "0") or "0"
        try:
            n_variants = int(n_var_str)
        except ValueError:
            n_variants = 0

        pgp_id = row.get("PGS Publication/Study (PGP) ID", "")
        pub = pubs.get(pgp_id, {})

        trait = efo_labels.split("|")[0].strip() if efo_labels else reported
        efo_id = efo_ids.split("|")[0].strip() if efo_ids else ""

        metadata[pgs_id] = {
            "trait": trait,
            "trait_efo": efo_id,
            "trait_reported": reported,
            "categories": [],
            "publication": pub.get("author", ""),
            "pub_year": pub.get("year", ""),
            "n_variants": n_variants,
        }

with open(output, "w") as f:
    json.dump(metadata, f, indent=2)

print(f"Wrote metadata for {len(metadata)} scores to {output}", file=sys.stderr)
