# Genome Query Skill

Query your personal genome (WGS VCF) for variants, genotypes, and annotations.

## Quick Start

```
/genome variant rs1805007      # Look up variant genotype + annotation
/genome region MC1R            # All variants in a gene
/genome region 16:89985000-89987000  # All variants in coordinates
/genome trait "red hair"       # GWAS Catalog trait lookup
/genome clinvar                # Scan for pathogenic variants
/genome hirisplex              # HIrisPlex-S 41-SNP panel
/genome lof                    # Loss-of-function variant scan
/genome impact MC1R            # All variants in gene with effect scores
/genome stats                  # Basic VCF statistics
/genome prs                    # Top/bottom polygenic risk scores
/genome prs intelligence       # Search PRS by trait name
/genome prs detail PGS003724   # Detailed view of single PRS
```

## Setup

### Environment Variable

Set `GENOME_DIR` to point to your genome data directory:

```bash
export GENOME_DIR=~/genome
```

### VCF Files

Place your VCF file(s) in `$GENOME_DIR` with these filenames:

| Build | Filename | Priority |
|-------|----------|----------|
| GRCh38 (annotated) | `grch38/genome.grch38.annotated.vcf.gz` | 1st |
| GRCh38 | `grch38/genome.grch38.vcf.gz` | 2nd |
| GRCh37 (rsID) | `genome.rsid.vcf.gz` | 3rd |
| GRCh37 (original) | `genome.vcf.gz` | 4th |

The skill auto-detects the best available VCF.

### PRS Results

If you've computed polygenic risk scores with pgsc_calc, place the results at:
`$GENOME_DIR/prs_results.json`

## Usage

### `/genome variant <rsID>`

Look up a single variant by rsID. Returns genotype from VCF + annotation from MyVariant.info.

```python
from genome_query import get_variant_by_rsid, annotate_variant

# Get genotype
result = get_variant_by_rsid("rs1805007")
# Returns: {'chrom': '16', 'pos': 89986117, 'id': 'rs1805007', 'ref': 'C', 'alt': 'T', 'gt': '0/0'}

# Get annotation
ann = annotate_variant("rs1805007")
# Returns ClinVar, gnomAD, CADD, SIFT, PolyPhen data
```

### `/genome region <gene_or_coords>`

Extract all variants in a gene or coordinate range.

```python
from genome_query import get_gene_variants, get_region

# By gene symbol (uses Ensembl REST API)
variants = get_gene_variants("MC1R")

# By coordinates
variants = get_region("16", 89985000, 89987000)
```

### `/genome trait "<trait>"`

Look up GWAS Catalog associations for a trait, then check your genotype for associated variants.

### `/genome clinvar`

Scan for pathogenic/likely pathogenic variants in your genome. Filters by ClinVar review status.

```python
from genome_query import scan_clinvar
results = scan_clinvar(min_review_stars=1)
```

### `/genome hirisplex`

Extract the HIrisPlex-S 41-SNP panel for hair/eye/skin colour prediction.

```python
from genome_query import extract_panel
results = extract_panel("hirisplex")
```

### `/genome lof`

Scan for putative loss-of-function variants:
- stop_gained, frameshift, splice donor/acceptor
- Missense with high CADD (>20) + deleterious SIFT + damaging PolyPhen
- Filtered by gnomAD AF <1% (rare variants more likely to be impactful)

```python
from genome_query import scan_lof
results = scan_lof(cadd_threshold=20, af_max=0.01)
```

### `/genome impact <gene>`

Show all variants in a gene with predicted effect scores (CADD, SIFT, PolyPhen, REVEL).

```python
from genome_query import get_gene_impact
results = get_gene_impact("BRCA1")
```

### `/genome stats`

Basic VCF statistics: variant count, Ti/Tv ratio, het/hom counts.

### `/genome prs`

Show top/bottom polygenic risk scores ranked by ancestry-adjusted percentile.

### `/genome prs <query>`

Search PRS results by trait name, EFO ID, PGS ID, or category.

### `/genome prs detail <PGS_ID>`

Detailed view of a single polygenic risk score.

## Data Sources

### Annotation Sources (via biothings-client → MyVariant.info)
- **ClinVar:** Clinical significance, review status, conditions
- **gnomAD:** Population allele frequencies (genome)
- **CADD:** Combined Annotation Dependent Depletion scores
- **SIFT:** Sorting Intolerant From Tolerant (missense prediction)
- **PolyPhen-2:** Polymorphism Phenotyping (missense prediction)
- **REVEL:** Rare Exome Variant Ensemble Learner
- **GWAS Catalog:** Trait associations
- **dbSNP:** rsID, allele info, MAF

### Gene Coordinates (via Ensembl REST API)
- **GRCh37:** `https://grch37.rest.ensembl.org/lookup/symbol/homo_sapiens/{gene}`
- **GRCh38:** `https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene}`

### Panels
- `panels/hirisplex.json` — 41 SNPs for hair/eye/skin colour (HIrisPlex-S)
- `panels/mc1r.json` — 8 MC1R red-hair variants

## Privacy Note

**Local-only operations (no external calls):**
- VCF queries via bcftools
- rsID annotation (bcftools + local dbSNP)
- All PRS queries — reads local JSON only

**External API calls (sends rsIDs/HGVS to MyVariant.info):**
- `annotate_variant()` — individual variant annotation
- `scan_lof()` — batch annotation of coding variants
- `get_gene_impact()` — batch annotation of gene variants
- Trait lookups via GWAS Catalog

MyVariant.info is a free, public bioinformatics API. No API key needed. No personal identifiers are sent — only variant IDs.

## Dependencies

- Python: `cyvcf2`, `biothings-client`
- System: `bcftools`, `tabix` (from htslib)
