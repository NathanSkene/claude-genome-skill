# Claude Genome Skill

A [Claude Code](https://docs.anthropic.com/en/docs/claude-code) skill for querying, annotating, and analysing your personal genome. Works with any whole genome sequencing (WGS) VCF file.

Built for the [OpenClaw Hackathon](https://github.com/anthropics/claude-code) — contributions welcome!

## Status

| Component | Status | Notes |
|-----------|--------|-------|
| **Query engine** (`genome_query.py`) | Tested | Used daily for personal genome queries |
| **SNP panels** (HIrisPlex, MC1R) | Tested | Verified against published genotype data |
| **SNP panels** (ACMG, pigmentation, oligodontia) | Tested | Used in personal analysis |
| **HPC pipeline scripts** | WIP | Developed for Imperial CX3 (PBS Pro). Paths are parameterised but scripts assume a specific HPC environment. You will likely need to adapt them. |
| **Triage pipeline** | WIP | vcfanno + slivar filtering works but the reference data download scripts need testing on fresh installs |
| **Report generator** | WIP | Generates HTML/Markdown reports. Template may need adjustment for different VCF formats |
| **PRS batch processing** | WIP | Requires pgsc_calc output. Batch retry logic tested on CX3 but not other schedulers |
| **Validation tools** | WIP | IGV-reports integration and variant checklist. Requires IGV-reports installation |

**"Tested"** = used repeatedly on real data. **"WIP"** = developed alongside a personal project, works in that context, but has not been systematically tested as a standalone package. Expect to read the scripts and adapt paths/configs for your setup.

## What it does

### Core (query engine)

- **Variant lookup** — Query any SNP by rsID and get your genotype + clinical annotation
- **Gene scanning** — Extract all variants in a gene with predicted effect scores (CADD, SIFT, PolyPhen)
- **ClinVar screening** — Scan your genome for pathogenic/likely pathogenic variants
- **Loss-of-function scan** — Find rare, high-impact coding variants
- **SNP panels** — Extract genotypes for predefined panels (HIrisPlex-S, MC1R, ACMG v3.3, pigmentation, oligodontia)
- **Polygenic risk scores** — View pre-computed PRS results with ancestry-adjusted percentiles
- **GWAS trait lookup** — Check your genotype for trait-associated variants

All queries run locally against your VCF. Annotation uses free public APIs (MyVariant.info, Ensembl REST) — no API keys needed.

### Extended pipeline (WIP)

- **HPC scripts** — End-to-end pipeline: BAM→FASTQ conversion, nf-core/raredisease variant calling, reference staging, container management
- **Standalone VEP annotation** — VEP v114 with CADD, LOFTEE, SpliceAI, AlphaMissense, and REVEL plugins (runs independently from the nf-core pipeline for easier upgrades)
- **Variant triage** — vcfanno annotation + slivar tiered filtering (artifact detection, allele frequency, clinical significance)
- **Clinical reports** — HTML + Markdown report generation from triage output
- **Variant validation** — Mini-BAM extraction for IGV, HTML report generation, interactive validation checklist

> **PRS computation** has moved to [vcf-to-prs](https://github.com/NathanSkene/vcf-to-prs). The PRS scripts in this repo (`prs_batch.sh`, `run_pgscalc.sh`, etc.) are kept for reference but are no longer maintained here.

## Quick Start

### 1. Get your genome

You need a VCF file from whole genome sequencing.

**Buy your own sequencing:**
- [Dante Labs](https://www.dantelabs.com/) (~$200, 30x WGS)
- [Nebula Genomics](https://nebula.org/) (~$250, 30x WGS)
- [23andMe](https://www.23andme.com/) (export VCF from settings — genotyping array, not WGS)

**Or test with public genomes (no purchase needed):**
- [Personal Genome Project (PGP)](https://www.personalgenomes.org/) — Hundreds of volunteers have made their WGS data fully public. Download VCFs from the [PGP Harvard data portal](https://my.pgp-hms.org/public_genetic_data). Try [hu661BA7 (GRCh38)](https://my.pgp-hms.org/profile/hu661BA7) — George Church's genome, one of the most studied public genomes.
- [1000 Genomes Project](https://www.internationalgenome.org/) — Population-scale WGS data. Individual sample VCFs available from the [IGSR data portal](https://www.internationalgenome.org/data-portal/sample). Try sample NA12878 (well-characterised reference genome).
- [openSNP](https://opensnp.org/) — Community-shared genotyping data (mostly array-based, not WGS, but works for variant lookups and panel extraction).

### 2. Generate your VCF (from FASTQ or BAM)

Most sequencing providers deliver raw data as FASTQ or BAM files, not VCFs. You need to run a variant calling pipeline.

**Recommended: [nf-core/raredisease](https://nf-co.re/raredisease)**

This is what we actually used. Raredisease handles the full workflow in one pipeline: alignment → variant calling (DeepVariant) → VEP annotation → clinical ranking. The triage and report scripts in this repo were built to consume its output.

```bash
nextflow run nf-core/raredisease \
    --input samplesheet.csv \
    --genome GRCh38 \
    --analysis_type wgs \
    -profile docker
```

**Alternative: [nf-core/sarek](https://nf-co.re/sarek)** — If you only need variant calling without the clinical annotation layer, sarek is lighter weight and uses GATK HaplotypeCaller.

See [claude-nextflow-skill](https://github.com/NathanSkene/claude-nextflow-skill) for automated samplesheet generation and pipeline orchestration.

**Running on HPC:** If you're on Imperial's CX3, see [claude-imperial-hpc-skill](https://github.com/NathanSkene/claude-imperial-hpc-skill) for PBS Pro configuration, or use the scripts in `hpc/` (WIP — these are configured for nf-core/raredisease on PBS Pro).

### 3. Install skill dependencies

```bash
# Python packages
pip3 install cyvcf2 biothings-client

# System tools (macOS)
brew install bcftools htslib

# System tools (Ubuntu/Debian)
sudo apt install bcftools tabix
```

### 4. Prepare your VCF

```bash
export GENOME_DIR=~/genome
mkdir -p "$GENOME_DIR"
# Copy your VCF to $GENOME_DIR/genome.vcf.gz

# Index it
tabix -p vcf "$GENOME_DIR/genome.vcf.gz"

# (Optional but recommended) Add rsID annotations from dbSNP
bcftools annotate -a dbsnp.vcf.gz -c ID genome.vcf.gz -Oz -o genome.rsid.vcf.gz
tabix -p vcf genome.rsid.vcf.gz
```

### 5. Install the skill

```bash
git clone https://github.com/NathanSkene/claude-genome-skill.git ~/.claude/skills/genome
export GENOME_DIR=~/genome  # add to your shell profile
```

### 6. Use it

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
├── README.md                    # This file
├── SKILL.md                     # Claude Code skill definition
├── METHODS.md                   # Reproducible pipeline documentation
├── genome-triggers.md           # Auto-detection rules
├── genome_query.py              # Core query engine (tested)
├── setup.sh                     # Bootstrap script
├── requirements.txt             # Python dependencies
│
├── panels/                      # SNP panels (tested)
│   ├── hirisplex.json           #   41-SNP pigmentation (HIrisPlex-S)
│   ├── mc1r.json                #   8 MC1R red hair variants
│   ├── acmg_sf_v3.3.json       #   81-gene ACMG secondary findings
│   ├── pigmentation.json        #   10-gene pigmentation panel (GRCh38)
│   └── oligodontia.json         #   8-gene oligodontia/ectodermal dysplasia
│
├── hpc/                         # HPC pipeline scripts (WIP — PBS Pro)
│   ├── setup_hpc.sh             #   Environment setup
│   ├── run_pipeline.sh          #   nf-core/raredisease launcher
│   ├── stage_references.sh      #   Reference genome staging
│   ├── download_fastq.sh        #   BAM/FASTQ download
│   ├── bam2fastq.sh             #   BAM→FASTQ conversion job
│   ├── run_pgscalc.sh           #   Single PRS batch
│   ├── prs_batch.sh             #   Multi-batch PRS orchestrator
│   ├── prs_retry_failed.sh      #   Retry failed PRS batches
│   ├── run_vep.sh                #   Standalone VEP v114 annotation
│   ├── monitor.sh               #   Job monitoring
│   ├── pull_containers.sh       #   Container download
│   ├── download_triage_refs.sh  #   Triage reference data
│   ├── download_gnomad_af.sh    #   gnomAD allele frequencies
│   ├── nextflow.config          #   PBS Pro + Apptainer config
│   ├── samplesheet.csv          #   nf-core input template
│   ├── prs_samplesheet.csv      #   pgsc_calc input template
│   ├── convert_metadata.py      #   Metadata format conversion
│   └── download_pgs_ids.py      #   PGS Catalog ID fetcher
│
├── triage/                      # Variant triage pipeline (WIP)
│   ├── run_triage.sh            #   vcfanno → slivar → tiering
│   ├── slivar_filters.js        #   Artifact detection + tiering
│   ├── vcfanno.toml             #   Annotation config
│   └── kosmos_check.sh          #   Known artifact comparison
│
├── report/                      # Clinical report generator (WIP)
│   ├── generate_report.py       #   HTML + Markdown report
│   └── report_template.html     #   Report template
│
├── prs/                         # PRS processing (WIP)
│   └── process_prs_results.py   #   pgsc_calc output → JSON
│
└── validation/                  # Variant validation (WIP)
    ├── extract_minibam.sh       #   Mini-BAM for IGV
    ├── run_igv_reports.sh       #   IGV-reports HTML generation
    └── validation_checklist.py  #   Interactive validation CLI
```

## Supported Genome Builds

The skill auto-detects your genome build:

| Build | VCF filename | Ensembl API |
|-------|-------------|-------------|
| GRCh38 (preferred) | `genome.grch38.vcf.gz` | rest.ensembl.org |
| GRCh37/hg19 | `genome.rsid.vcf.gz` or `genome.vcf.gz` | grch37.rest.ensembl.org |

## Extended Pipeline Overview (WIP)

The full analysis pipeline runs in this order:

```
FASTQ/BAM → nf-core/raredisease → VCF → triage → report
                                    ↓
                              pgsc_calc → PRS results → genome_query.py
```

### HPC Pipeline (`hpc/`)

Scripts assume PBS Pro on Imperial CX3 but can be adapted for SLURM or other schedulers. All scripts use environment variables for paths:

```bash
export HPC_USER=your_username
export SAMPLE_ID=your_sample
export GENOME_DIR=/path/to/your/project
```

### Standalone VEP Annotation (`hpc/run_vep.sh`)

VEP annotation runs as a separate PBS job rather than inside nf-core/raredisease. This decouples the annotation version from the pipeline container, making upgrades straightforward.

**Current setup:** VEP v114 container with plugins:
- **CADD v1.7** — Combined deleteriousness score
- **LOFTEE** — Loss-of-function classification (HC/LC)
- **SpliceAI** — Splice disruption prediction
- **AlphaMissense** — Protein structure-based pathogenicity
- **REVEL** — Ensemble missense predictor

VEP v114 cache also includes MaveDB (7.7M functionally-assessed variants), PrimateAI-3D, updated LOEUF v4.1 constraint scores, and AllOfUs population frequencies via dbNSFP v5.0c.

```bash
# Stage all references (VEP cache, CADD, SpliceAI, AlphaMissense, REVEL, LOFTEE)
bash hpc/stage_references.sh

# Run VEP on DeepVariant output
qsub hpc/run_vep.sh /path/to/deepvariant.vcf.gz
```

### Triage Pipeline (`triage/`)

Three-stage variant filtering:
1. **vcfanno** — Annotate with segmental duplications, low-complexity regions, gnomAD AF
2. **slivar** — Artifact detection + tiered filtering (pathogenic, protein-altering, rare)
3. **Clinical tiering** — Tier 1 (pathogenic), Tier 2 (likely pathogenic), Tier 3 (VUS)

Requires: slivar, vcfanno, bcftools, and reference data (download with `hpc/download_triage_refs.sh`).

### Report Generator (`report/`)

Produces a clinical-style HTML report from triage output. Includes variant tables, gene summaries, and ACMG classification.

### PRS Batch System → [vcf-to-prs](https://github.com/NathanSkene/vcf-to-prs)

PRS computation has been moved to a dedicated repository. See [vcf-to-prs](https://github.com/NathanSkene/vcf-to-prs) for pgsc_calc orchestration, batch processing, and result parsing.

## Polygenic Risk Scores

PRS requires pre-computation using [pgsc_calc](https://pgsc-calc.readthedocs.io/):

```bash
nextflow run pgscatalog/pgsc_calc \
    --input samplesheet.csv \
    --pgs_id PGS003724,PGS003725 \
    --target_build GRCh37 \
    --run_ancestry \
    -profile singularity
```

The skill reads `$GENOME_DIR/prs_results.json`.

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

The WIP components need the most help. If you test any of them on your setup, please open an issue describing what worked and what needed changing.

Ideas for improvement:
- [ ] SLURM equivalents for the PBS Pro HPC scripts
- [ ] Docker/Singularity containerisation of the triage pipeline
- [ ] More SNP panels (pharmacogenomics, ancestry informative markers)
- [ ] GRCh38 liftover utility
- [ ] Integration with [ClinGen](https://clinicalgenome.org/) actionability
- [ ] Carrier status screening panels
- [ ] Web UI for results visualisation

## Related

- [vcf-to-prs](https://github.com/NathanSkene/vcf-to-prs) — Polygenic risk score computation pipeline (pgsc_calc orchestration, batch processing)
- [claude-imperial-hpc-skill](https://github.com/NathanSkene/claude-imperial-hpc-skill) — PBS Pro reference for Imperial's CX3 cluster

## License

MIT
