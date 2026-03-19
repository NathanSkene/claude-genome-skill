#!/usr/bin/env python3
"""
Process pgsc_calc output into a queryable JSON file with z-scores.

Reads concatenated scores from pgsc_calc, enriches with trait metadata,
downloads scoring files to compute theoretical population distributions,
and outputs prs_results.json with z-scores and percentiles.

Z-score method: For each PGS, downloads the scoring file from PGS Catalog FTP
to get effect weights and allele frequencies. Computes theoretical population
mean and SD under Hardy-Weinberg equilibrium for a diploid genome:
  mean = sum(2 * weight_i * af_i)
  sd   = sqrt(sum(2 * weight_i^2 * af_i * (1 - af_i)))
When allele frequencies are not in the scoring file, uses 0.5 (uninformative
prior) which gives approximate but still informative z-scores.

Caveats:
  - Population stats assume EUR ancestry (most PGS are EUR-calibrated)
  - AF=0.5 fallback is less accurate for variants with extreme MAF
  - Proper ancestry adjustment (PCA-based via pgsc_calc --run_ancestry) is
    more rigorous but requires internet + reference panels on compute nodes

Usage:
    python3 process_prs_results.py scores.txt.gz --metadata pgs_metadata.json
    python3 process_prs_results.py scores.txt.gz --no-zscore  # skip scoring file download
"""

import argparse
import csv
import gzip
import io
import json
import math
import os
import sys
import time
import urllib.request
import urllib.error
from datetime import date

PGS_API = "https://www.pgscatalog.org/rest"
PGS_FTP = "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores"
DEFAULT_AF = 0.5  # Fallback when scoring file lacks allele frequencies


def fetch_json(url: str, retries: int = 3) -> dict | None:
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
            print(f"Warning: Failed to fetch {url}: {e}", file=sys.stderr)
            return None


def fetch_bytes(url: str, retries: int = 3) -> bytes | None:
    """Fetch raw bytes from URL with retry logic."""
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url)
            with urllib.request.urlopen(req, timeout=60) as resp:
                return resp.read()
        except (urllib.error.URLError, urllib.error.HTTPError) as e:
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
                continue
            print(f"Warning: Failed to fetch {url}: {e}", file=sys.stderr)
            return None


def load_scores(scores_path: str) -> list[dict]:
    """Load pgsc_calc aggregated scores file."""
    rows = []
    opener = gzip.open if scores_path.endswith(".gz") else open

    with opener(scores_path, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)

    return rows


def load_metadata_file(metadata_path: str) -> dict:
    """Load pre-downloaded PGS metadata JSON."""
    with open(metadata_path) as f:
        return json.load(f)


def fetch_metadata_for_ids(pgs_ids: list[str]) -> dict:
    """Fetch trait metadata from PGS Catalog API for a list of IDs."""
    metadata = {}
    total = len(pgs_ids)

    for i, pgs_id in enumerate(pgs_ids):
        if (i + 1) % 50 == 0:
            print(f"  Fetching metadata: {i+1}/{total}", file=sys.stderr)

        data = fetch_json(f"{PGS_API}/score/{pgs_id}?format=json")
        if data is None:
            metadata[pgs_id] = {
                "trait": "Unknown",
                "trait_efo": "",
                "trait_reported": "",
                "categories": [],
                "publication": "",
                "pub_year": "",
                "n_variants": 0,
            }
            continue

        traits = data.get("trait_efo", [])
        trait_name = traits[0]["label"] if traits else data.get("trait_reported", "Unknown")
        trait_id = traits[0]["id"] if traits else ""
        categories = []
        for t in traits:
            categories.extend(t.get("trait_categories", []))

        pub = data.get("publication", {}) or {}
        metadata[pgs_id] = {
            "trait": trait_name,
            "trait_efo": trait_id,
            "trait_reported": data.get("trait_reported", ""),
            "categories": list(set(categories)),
            "publication": pub.get("firstauthor", ""),
            "pub_year": pub.get("date_publication", "")[:4] if pub.get("date_publication") else "",
            "n_variants": data.get("variants_number", 0),
        }

        time.sleep(0.2)  # Rate limiting

    return metadata


def parse_scoring_file(data: bytes) -> dict:
    """Parse a PGS Catalog scoring file (gzipped) into weights and allele frequencies.

    Returns dict with:
        weights: list of float effect weights
        afs: list of float allele frequencies (DEFAULT_AF if not in file)
        n_variants: number of data rows
        has_af: whether the file contained allele frequency data
    """
    try:
        text = gzip.decompress(data).decode("utf-8")
    except Exception:
        return {"weights": [], "afs": [], "n_variants": 0, "has_af": False}

    lines = text.split("\n")

    # Find header line (first non-comment, non-empty line)
    header_line = None
    data_start = 0
    for i, line in enumerate(lines):
        if line.strip() and not line.startswith("#"):
            header_line = line
            data_start = i + 1
            break

    if header_line is None:
        return {"weights": [], "afs": [], "n_variants": 0, "has_af": False}

    cols = header_line.strip().split("\t")
    col_idx = {c.lower(): i for i, c in enumerate(cols)}

    weight_col = col_idx.get("effect_weight")
    if weight_col is None:
        return {"weights": [], "afs": [], "n_variants": 0, "has_af": False}

    # Check for allele frequency columns (various naming conventions)
    af_col = None
    for af_name in ["allelefrequency_effect", "effect_allele_frequency",
                     "allele_frequency", "eaf", "af"]:
        if af_name in col_idx:
            af_col = col_idx[af_name]
            break

    weights = []
    afs = []
    for line in lines[data_start:]:
        if not line.strip():
            continue
        fields = line.strip().split("\t")
        if len(fields) <= weight_col:
            continue

        try:
            w = float(fields[weight_col])
        except (ValueError, IndexError):
            continue

        af = DEFAULT_AF
        if af_col is not None and len(fields) > af_col:
            try:
                af = float(fields[af_col])
                if af <= 0 or af >= 1:
                    af = DEFAULT_AF
            except (ValueError, IndexError):
                af = DEFAULT_AF

        weights.append(w)
        afs.append(af)

    return {
        "weights": weights,
        "afs": afs,
        "n_variants": len(weights),
        "has_af": af_col is not None,
    }


def compute_population_stats(weights: list[float], afs: list[float]) -> dict:
    """Compute theoretical population mean and SD for a PGS under HWE (diploid).

    For diploid genome with allele frequencies p_i and effect weights w_i:
        E[PRS] = sum(2 * w_i * p_i)
        Var[PRS] = sum(2 * w_i^2 * p_i * (1 - p_i))
        SD[PRS] = sqrt(Var[PRS])
    """
    if not weights:
        return {"mean": None, "sd": None}

    mean = sum(2 * w * p for w, p in zip(weights, afs))
    variance = sum(2 * w * w * p * (1 - p) for w, p in zip(weights, afs))
    sd = math.sqrt(variance) if variance > 0 else None

    return {"mean": mean, "sd": sd}


def risk_category(percentile: float) -> str:
    """Assign risk category based on percentile."""
    if percentile < 20:
        return "Low"
    elif percentile <= 80:
        return "Average"
    elif percentile <= 95:
        return "Elevated"
    else:
        return "High"


def normal_cdf(z: float) -> float:
    """Approximate the standard normal CDF using the error function.

    Uses math.erf which is available in Python 3.2+, avoiding scipy dependency.
    """
    return 0.5 * (1 + math.erf(z / math.sqrt(2)))


def fetch_population_stats(pgs_ids: list[str], build: str,
                           cache_path: str | None = None) -> dict:
    """Download scoring files and compute population distribution stats.

    Args:
        pgs_ids: list of base PGS IDs (e.g. PGS000033, not PGS000033_hmPOS_GRCh37)
        build: genome build (GRCh37 or GRCh38)
        cache_path: optional path to cache scoring file stats as JSON

    Returns:
        dict mapping PGS ID to {mean, sd, n_variants, has_af}
    """
    # Load cache if it exists
    cached = {}
    if cache_path and os.path.exists(cache_path):
        with open(cache_path) as f:
            cached = json.load(f)
        print(f"  Loaded {len(cached)} cached population stats", file=sys.stderr)

    stats = {}
    to_fetch = [p for p in pgs_ids if p not in cached]
    total = len(to_fetch)

    if cached:
        for p in pgs_ids:
            if p in cached:
                stats[p] = cached[p]

    if to_fetch:
        print(f"  Downloading {total} scoring files from PGS Catalog FTP...",
              file=sys.stderr)

    for i, pgs_id in enumerate(to_fetch):
        if (i + 1) % 10 == 0 or (i + 1) == total:
            print(f"  Scoring files: {i+1}/{total}", file=sys.stderr)

        url = f"{PGS_FTP}/{pgs_id}/ScoringFiles/Harmonized/{pgs_id}_hmPOS_{build}.txt.gz"
        data = fetch_bytes(url)

        if data is None:
            # Try unharmonized
            url = f"{PGS_FTP}/{pgs_id}/ScoringFiles/{pgs_id}.txt.gz"
            data = fetch_bytes(url)

        if data is None:
            stats[pgs_id] = {"mean": None, "sd": None, "n_variants": 0,
                             "has_af": False}
            continue

        parsed = parse_scoring_file(data)
        pop_stats = compute_population_stats(parsed["weights"], parsed["afs"])

        stats[pgs_id] = {
            "mean": pop_stats["mean"],
            "sd": pop_stats["sd"],
            "n_variants": parsed["n_variants"],
            "has_af": parsed["has_af"],
        }

        time.sleep(0.1)  # Rate limiting

    # Save cache
    if cache_path:
        all_stats = {**cached, **{p: stats[p] for p in to_fetch}}
        with open(cache_path, "w") as f:
            json.dump(all_stats, f, indent=2)
        print(f"  Cached population stats to {cache_path}", file=sys.stderr)

    return stats


def process_scores(scores: list[dict], metadata: dict,
                   pop_stats: dict | None = None) -> list[dict]:
    """Process raw scores into enriched result entries with z-scores."""
    results = []

    for row in scores:
        pgs_id = row.get("PGS") or row.get("PGS_ID") or row.get("pgs_id", "")
        if not pgs_id:
            continue

        # Strip pgsc_calc suffix (e.g., PGS000033_hmPOS_GRCh37 -> PGS000033)
        base_id = pgs_id.split("_")[0] if "_" in pgs_id else pgs_id
        meta = metadata.get(base_id, metadata.get(pgs_id, {}))

        raw_sum = _safe_float(row.get("SUM"))

        entry = {
            "pgs_id": base_id,
            "trait": meta.get("trait", "Unknown"),
            "trait_efo": meta.get("trait_efo", ""),
            "trait_reported": meta.get("trait_reported", ""),
            "category": ", ".join(meta.get("categories", [])) if meta.get("categories") else "Uncategorized",
            "sum": raw_sum,
        }

        # Variant counts from scoring file stats (accurate) or metadata (fallback)
        n_total = meta.get("n_variants", 0)
        if pop_stats and base_id in pop_stats:
            ps = pop_stats[base_id]
            n_in_file = ps.get("n_variants", 0)
            if n_in_file > 0:
                n_total = n_in_file
        entry["n_variants_total"] = n_total

        # Z-scores from theoretical population distribution
        if pop_stats and base_id in pop_stats:
            ps = pop_stats[base_id]
            pop_mean = ps.get("mean")
            pop_sd = ps.get("sd")
            has_af = ps.get("has_af", False)

            if raw_sum is not None and pop_mean is not None and pop_sd and pop_sd > 0:
                z = (raw_sum - pop_mean) / pop_sd
                pct = normal_cdf(z) * 100
                entry["z_score"] = round(z, 4)
                entry["percentile"] = round(pct, 2)
                entry["risk_category"] = risk_category(pct)
                entry["has_population_af"] = has_af
            else:
                entry["z_score"] = None
                entry["percentile"] = None
                entry["risk_category"] = None
                entry["has_population_af"] = has_af
        else:
            # Ancestry-adjusted columns from pgsc_calc --run_ancestry (if present)
            if "percentile_MostSimilarPop" in row:
                entry["percentile"] = _safe_float(row.get("percentile_MostSimilarPop"))
            if "Z_MostSimilarPop" in row:
                entry["z_score"] = _safe_float(row.get("Z_MostSimilarPop"))
            if "MostSimilarPop" in row:
                entry["ancestry_pop"] = row.get("MostSimilarPop", "")

        # Publication info
        entry["publication"] = meta.get("publication", "")
        entry["pub_year"] = meta.get("pub_year", "")

        results.append(entry)

    return results


def _safe_float(val) -> float | None:
    """Convert to float or None."""
    if val is None or val == "" or val == "NA":
        return None
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


def _safe_int(val) -> int | None:
    """Convert to int or None."""
    if val is None or val == "" or val == "NA":
        return None
    try:
        return int(float(val))
    except (ValueError, TypeError):
        return None


def main():
    parser = argparse.ArgumentParser(description="Process pgsc_calc PRS output")
    parser.add_argument("scores",
                        help="Path to pgsc_calc aggregated scores file (.txt.gz or .txt)")
    parser.add_argument("--metadata", "-m",
                        help="Pre-downloaded PGS metadata JSON (from download_pgs_ids.py)")
    parser.add_argument("--output", "-o", default="prs_results.json",
                        help="Output JSON file")
    parser.add_argument("--build", default="GRCh37",
                        help="Genome build used (default: GRCh37)")
    parser.add_argument("--sample", default="YOUR_SAMPLE",
                        help="Sample ID")
    parser.add_argument("--no-zscore", action="store_true",
                        help="Skip scoring file download and z-score computation")
    parser.add_argument("--stats-cache", default="pgs_population_stats.json",
                        help="Cache file for population stats (default: pgs_population_stats.json)")
    args = parser.parse_args()

    # Load scores
    print(f"Loading scores from {args.scores}...", file=sys.stderr)
    scores = load_scores(args.scores)
    print(f"  Loaded {len(scores)} score rows", file=sys.stderr)

    # Get unique base PGS IDs
    pgs_ids_raw = list(set(
        row.get("PGS") or row.get("PGS_ID") or row.get("pgs_id", "")
        for row in scores
    ))
    pgs_ids_raw = [p for p in pgs_ids_raw if p]
    base_ids = sorted(set(p.split("_")[0] if "_" in p else p for p in pgs_ids_raw))
    print(f"  Unique PGS IDs: {len(base_ids)}", file=sys.stderr)

    # Load or fetch metadata
    if args.metadata and os.path.exists(args.metadata):
        print(f"Loading metadata from {args.metadata}...", file=sys.stderr)
        metadata = load_metadata_file(args.metadata)
        # Fetch any missing IDs
        missing = [b for b in base_ids if b not in metadata]
        if missing:
            print(f"  Fetching metadata for {len(missing)} missing IDs...",
                  file=sys.stderr)
            extra = fetch_metadata_for_ids(missing)
            metadata.update(extra)
    else:
        print(f"Fetching metadata for {len(base_ids)} PGS IDs from API...",
              file=sys.stderr)
        metadata = fetch_metadata_for_ids(base_ids)

    # Fetch population distribution stats from scoring files
    pop_stats = None
    if not args.no_zscore:
        print("Computing population distribution stats...", file=sys.stderr)
        pop_stats = fetch_population_stats(
            base_ids, args.build, cache_path=args.stats_cache)
        with_af = sum(1 for s in pop_stats.values() if s.get("has_af"))
        with_sd = sum(1 for s in pop_stats.values()
                      if s.get("sd") is not None and s["sd"] > 0)
        print(f"  {with_sd}/{len(pop_stats)} have computable z-scores "
              f"({with_af} with published AFs, "
              f"{with_sd - with_af} using AF=0.5 approximation)",
              file=sys.stderr)

    # Process
    print("Processing scores...", file=sys.stderr)
    results = process_scores(scores, metadata, pop_stats)

    # Sort by percentile (highest risk first), then by absolute z-score
    results.sort(
        key=lambda x: (
            x.get("percentile") is not None,  # scored first
            x.get("percentile") or 0,
        ),
        reverse=True,
    )

    # Build output
    n_with_z = sum(1 for r in results if r.get("z_score") is not None)
    output = {
        "sample": args.sample,
        "build": args.build,
        "n_scores": len(results),
        "n_with_zscore": n_with_z,
        "zscore_method": "theoretical_population" if not args.no_zscore else None,
        "zscore_note": (
            "Z-scores computed from theoretical population distribution under "
            "Hardy-Weinberg equilibrium. Scores without published allele frequencies "
            "use AF=0.5 approximation (has_population_af=false). "
            "Percentiles assume European ancestry."
        ) if not args.no_zscore else None,
        "date_computed": str(date.today()),
        "scores": results,
    }

    # Write
    with open(args.output, "w") as f:
        json.dump(output, f, indent=2)

    print(f"\nWrote {len(results)} scores to {args.output} "
          f"({n_with_z} with z-scores)", file=sys.stderr)

    # Summary stats
    with_z = [r for r in results if r.get("z_score") is not None]
    if with_z:
        print("\nRisk summary:", file=sys.stderr)
        cats = {}
        for r in with_z:
            cat = r.get("risk_category", "Unknown")
            cats[cat] = cats.get(cat, 0) + 1
        for cat in ["High", "Elevated", "Average", "Low"]:
            if cat in cats:
                print(f"  {cat}: {cats[cat]}", file=sys.stderr)

        # Top 5 highest risk
        high_risk = sorted(with_z, key=lambda x: x["percentile"], reverse=True)
        print("\nTop 5 highest percentile:", file=sys.stderr)
        for r in high_risk[:5]:
            af_note = "" if r.get("has_population_af") else " (AF≈0.5)"
            print(f"  {r['pgs_id']} ({r['trait']}): "
                  f"z={r['z_score']:+.2f}, "
                  f"{r['percentile']:.1f}th pctl, "
                  f"{r['risk_category']}{af_note}",
                  file=sys.stderr)

        # Top 5 lowest risk
        print("\nTop 5 lowest percentile:", file=sys.stderr)
        low_risk = sorted(with_z, key=lambda x: x["percentile"])
        for r in low_risk[:5]:
            af_note = "" if r.get("has_population_af") else " (AF≈0.5)"
            print(f"  {r['pgs_id']} ({r['trait']}): "
                  f"z={r['z_score']:+.2f}, "
                  f"{r['percentile']:.1f}th pctl, "
                  f"{r['risk_category']}{af_note}",
                  file=sys.stderr)


if __name__ == "__main__":
    main()
