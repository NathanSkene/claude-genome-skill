# Claude Genome Skill

A [Claude Code](https://docs.anthropic.com/en/docs/claude-code) skill for querying and annotating your personal genome. Works with any whole genome sequencing (WGS) VCF file.

Built for the [OpenClaw Hackathon](https://github.com/anthropics/claude-code) — contributions welcome!

## What it does

- **Variant lookup** — Query any SNP by rsID and get your genotype + clinical annotation
- **Gene scanning** — Extract all variants in a gene with predicted effect scores (CADD, SIFT, PolyPhen)
- **ClinVar screening** — Scan your genome for pathogenic/likely pathogenic variants
- **Loss-of-function scan** — Find rare, high-impact coding variants
- **SNP panels** — Extract genotypes for predefined panels (HIrisPlex-S pigmentation, MC1R red hair)
- **Polygenic risk scores** — View pre-computed PRS results with ancestry-adjusted percentiles
- **GWAS trait lookup** — Check your genotype for trait-associated variants

All queries run locally against your VCF. Annotation uses free public APIs (MyVariant.info, Ensembl REST) — no API keys needed.

## Quick Start

### 1. Get your genome

You need a VCF file from whole genome sequencing. Providers:
- [Dante Labs](https://www.dantelabs.com/) (~$200, 30x WGS)
- [Nebula Genomics](https://nebula.org/) (~$250, 30x WGS)
- [23andMe](https://www.23andme.com/) (export VCF from settings — genotyping array, not WGS)

### 2. Install dependencies

```bash
# Python packages
pip3 install cyvcf2 biothings-client

# System tools (macOS)
brew install bcftools htslib

# System tools (Ubuntu/Debian)
sudo apt install bcftools tabix
```

### 3. Prepare your VCF

```bash
# Set your genome directory
export GENOME_DIR=~/genome

mkdir -p "$GENOME_DIR"
# Copy your VCF to $GENOME_DIR/genome.vcf.gz

# Index it
tabix -p vcf "$GENOME_DIR/genome.vcf.gz"

# (Optional but recommended) Add rsID annotations from dbSNP
# Download dbSNP VCF for your build (GRCh37 or GRCh38)
# Then annotate:
bcftools annotate -a dbsnp.vcf.gz -c ID genome.vcf.gz -Oz -o genome.rsid.vcf.gz
tabix -p vcf genome.rsid.vcf.gz
```

### 4. Install the skill

```bash
# Clone to your Claude Code skills directory
git clone https://github.com/NathanSkene/claude-genome-skill.git ~/.claude/skills/genome

# Set environment variable (add to your shell profile)
export GENOME_DIR=~/genome
```

### 5. Use it

In Claude Code:
```
/genome variant rs1805007      # MC1R red hair variant
/genome region BRCA1           # All variants in BRCA1
/genome hirisplex              # Hair/eye/skin colour prediction
/genome clinvar                # Pathogenic variant scan
/genome lof                    # Loss-of-function variants
/genome prs                    # Polygenic risk scores (if computed)
```

## File Structure

```
claude-genome-skill/
├── README.md              # This file
├── SKILL.md               # Claude Code skill definition
├── genome-triggers.md     # Auto-detection rules (copy to ~/.claude/rules/)
├── genome_query.py        # Core query engine
├── panels/
│   ├── hirisplex.json     # 41-SNP pigmentation panel
│   └── mc1r.json          # 8 MC1R red hair variants
├── setup.sh               # Bootstrap script
└── requirements.txt       # Python dependencies
```

## Supported Genome Builds

The skill auto-detects your genome build:

| Build | VCF filename | Ensembl API |
|-------|-------------|-------------|
| GRCh38 (preferred) | `genome.grch38.vcf.gz` | rest.ensembl.org |
| GRCh37/hg19 | `genome.rsid.vcf.gz` or `genome.vcf.gz` | grch37.rest.ensembl.org |

Place your VCF in `$GENOME_DIR` with the appropriate filename.

## Polygenic Risk Scores

PRS requires pre-computation using [pgsc_calc](https://pgsc-calc.readthedocs.io/):

```bash
# Run on HPC (needs 16-128GB RAM depending on number of scores)
nextflow run pgscatalog/pgsc_calc \
    --input samplesheet.csv \
    --pgs_id PGS003724,PGS003725 \
    --target_build GRCh37 \
    --run_ancestry \
    -profile singularity
```

The skill reads `$GENOME_DIR/prs_results.json`. If you're on Imperial's CX3 cluster, see [claude-imperial-hpc-skill](https://github.com/NathanSkene/claude-imperial-hpc-skill) for PBS Pro configuration.

## Adding Custom Panels

Create a JSON file in `panels/`:

```json
{
  "name": "My Custom Panel",
  "description": "Description of the panel",
  "reference": "GRCh37",
  "snps": [
    {
      "rsid": "rs12345",
      "gene": "GENE",
      "chrom": "1",
      "pos": 12345678,
      "ref": "A",
      "alt": "G",
      "trait": "trait description"
    }
  ]
}
```

Then: `/genome panel my_custom_panel`

## Privacy

**Local-only operations** (no data leaves your machine):
- All VCF queries (variant lookup, region extraction, panels, PRS)

**External API calls** (sends rsIDs to public databases):
- `annotate_variant()` → MyVariant.info
- `scan_clinvar()` → MyVariant.info
- `scan_lof()` → MyVariant.info
- `get_gene_impact()` → MyVariant.info + Ensembl REST
- Gene coordinate lookup → Ensembl REST

Only variant IDs are sent (e.g., "rs1805007"). No personal identifiers.

## Contributing

Ideas for improvement:
- [ ] More SNP panels (pharmacogenomics, ancestry informative markers, ACMG73)
- [ ] GRCh38 liftover utility
- [ ] Integration with [ClinGen](https://clinicalgenome.org/) actionability
- [ ] Carrier status screening panels
- [ ] Batch VCF comparison (family/trio analysis)
- [ ] Web UI for results visualization

## Related

- [claude-imperial-hpc-skill](https://github.com/NathanSkene/claude-imperial-hpc-skill) — PBS Pro reference for Imperial's CX3 cluster

## License

MIT
