#!/usr/bin/env python3
"""
Download all PGS Catalog score IDs via REST API.

Outputs a text file with one PGS ID per line, optionally filtered by
trait category or EFO term.

Usage:
    python3 download_pgs_ids.py                          # All scores
    python3 download_pgs_ids.py --category "Body measurement"
    python3 download_pgs_ids.py --efo EFO_0004337        # Intelligence
    python3 download_pgs_ids.py --output pgs_ids.txt
"""

import argparse
import json
import sys
import time
import urllib.request
import urllib.error

PGS_API = "https://www.pgscatalog.org/rest"


def fetch_json(url: str, retries: int = 3) -> dict:
    """Fetch JSON from URL with retry logic."""
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                return json.loads(resp.read().decode())
        except (urllib.error.URLError, urllib.error.HTTPError) as e:
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
                continue
            raise RuntimeError(f"Failed to fetch {url}: {e}")


def get_all_score_ids():
    """Fetch all PGS IDs from the catalog, paginating through results."""
    ids = []
    url = f"{PGS_API}/score/all?format=json&limit=250"

    while url:
        print(f"Fetching: {url}", file=sys.stderr)
        data = fetch_json(url)
        for score in data.get("results", []):
            ids.append(score["id"])
        url = data.get("next")
        if url:
            time.sleep(0.5)  # Rate limiting

    return ids


def get_score_ids_by_trait(efo_id):
    """Fetch PGS IDs for a specific EFO trait."""
    ids = []
    url = f"{PGS_API}/score/search?trait_id={efo_id}&format=json&limit=250"

    while url:
        print(f"Fetching: {url}", file=sys.stderr)
        data = fetch_json(url)
        for score in data.get("results", []):
            ids.append(score["id"])
        url = data.get("next")
        if url:
            time.sleep(0.5)

    return ids


def get_score_ids_by_category(category):
    """Fetch PGS IDs filtered by trait category name (substring match)."""
    all_ids = []
    url = f"{PGS_API}/score/all?format=json&limit=250"

    while url:
        print(f"Fetching: {url}", file=sys.stderr)
        data = fetch_json(url)
        for score in data.get("results", []):
            trait_efo = score.get("trait_efo", [])
            for trait in trait_efo:
                cat = trait.get("trait_categories", [])
                if any(category.lower() in c.lower() for c in cat):
                    all_ids.append(score["id"])
                    break
        url = data.get("next")
        if url:
            time.sleep(0.5)

    return all_ids


def get_trait_metadata(pgs_ids):
    """Fetch trait metadata for a list of PGS IDs (for enrichment)."""
    metadata = {}
    total = len(pgs_ids)
    for i, pgs_id in enumerate(pgs_ids):
        if (i + 1) % 100 == 0:
            print(f"  Metadata: {i+1}/{total}", file=sys.stderr)
        url = f"{PGS_API}/score/{pgs_id}?format=json"
        try:
            data = fetch_json(url)
            traits = data.get("trait_efo", [])
            trait_name = traits[0]["label"] if traits else data.get("trait_reported", "Unknown")
            trait_id = traits[0]["id"] if traits else ""
            categories = []
            for t in traits:
                categories.extend(t.get("trait_categories", []))

            metadata[pgs_id] = {
                "trait": trait_name,
                "trait_efo": trait_id,
                "trait_reported": data.get("trait_reported", ""),
                "categories": list(set(categories)),
                "publication": data.get("publication", {}).get("firstauthor", ""),
                "pub_year": data.get("publication", {}).get("date_publication", "")[:4] if data.get("publication", {}).get("date_publication") else "",
                "n_variants": data.get("variants_number", 0),
            }
            time.sleep(1.0)
        except Exception as e:
            print(f"Warning: Failed to fetch metadata for {pgs_id}: {e}", file=sys.stderr)
            metadata[pgs_id] = {"trait": "Unknown", "trait_efo": "", "categories": []}
            time.sleep(5.0)  # Back off on errors

    return metadata


def main():
    parser = argparse.ArgumentParser(description="Download PGS Catalog score IDs")
    parser.add_argument("--efo", help="Filter by EFO trait ID (e.g., EFO_0004337)")
    parser.add_argument("--category", help="Filter by trait category name")
    parser.add_argument("--output", "-o", default="pgs_ids.txt", help="Output file (default: pgs_ids.txt)")
    parser.add_argument("--metadata", action="store_true", help="Also fetch trait metadata as JSON")
    parser.add_argument("--metadata-output", default="pgs_metadata.json", help="Metadata output file")
    args = parser.parse_args()

    if args.efo:
        print(f"Fetching scores for trait: {args.efo}", file=sys.stderr)
        ids = get_score_ids_by_trait(args.efo)
    elif args.category:
        print(f"Fetching scores for category: {args.category}", file=sys.stderr)
        ids = get_score_ids_by_category(args.category)
    else:
        print("Fetching all PGS Catalog score IDs...", file=sys.stderr)
        ids = get_all_score_ids()

    # Sort IDs numerically
    ids.sort(key=lambda x: int(x.replace("PGS", "")))

    with open(args.output, "w") as f:
        for pgs_id in ids:
            f.write(pgs_id + "\n")

    print(f"Wrote {len(ids)} PGS IDs to {args.output}", file=sys.stderr)

    if args.metadata:
        print(f"Fetching trait metadata for {len(ids)} scores...", file=sys.stderr)
        metadata = get_trait_metadata(ids)
        with open(args.metadata_output, "w") as f:
            json.dump(metadata, f, indent=2)
        print(f"Wrote metadata to {args.metadata_output}", file=sys.stderr)


if __name__ == "__main__":
    main()
