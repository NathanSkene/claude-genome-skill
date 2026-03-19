# Genome Analysis Pipeline — Methods

Reproducible documentation for the genome analysis pipeline setup.

**Date:** 2026-02-22
**Author:** YOUR_NAME (with Claude Code)
**Sample:** YOUR_SAMPLE (Dante Labs WGS)

---

## 1. Source Data

- **VCF file:** `MY GENOME FROM DANTE LABS - COMBINED VCF - YOUR_SAMPLE.filtered.vcf.gz`
- **Original location:** OneDrive cloud storage
- **Reference genome:** GRCh37/hg19 (verified from VCF header `##reference=file:///references/hs37d5/hs37d5.fa`)
- **Chromosome naming:** Plain numbers (1, 2, ... — no "chr" prefix)
- **Variant count:** ~5.2M

## 2. Local Copy + Index

Copied VCF to local directory for performance (avoids OneDrive latency):

```bash
cd /path/to/GenomeAnalysis/
cp "/path/to/OneDrive/MY GENOME FROM DANTE LABS - COMBINED VCF - YOUR_SAMPLE.filtered.vcf.gz" genome.vcf.gz
tabix -p vcf genome.vcf.gz
```

**Tools:** tabix (from htslib, installed via Homebrew)

## 3. rsID Annotation with dbSNP

### 3a. Download full dbSNP b151

Using the full dbSNP VCF (not common_all) to avoid missing rare/ClinVar variants:

```bash
cd /path/to/GenomeAnalysis/
curl -O https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz
curl -O https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz.tbi
```

- **File:** `00-All.vcf.gz` (~15GB)
- **dbSNP build:** b151
- **Assembly:** GRCh37p13
- **Chromosome naming:** Plain numbers (matches VCF)

### 3b. Annotate VCF with rsIDs

```bash
bcftools annotate -a 00-All.vcf.gz -c ID -Oz -o genome.rsid.vcf.gz genome.vcf.gz
tabix -p vcf genome.rsid.vcf.gz
```

### 3c. Verify annotation

```bash
# Count rsID-annotated variants
bcftools query -f '%ID\n' genome.rsid.vcf.gz | grep -c "^rs"
# Expected: ~4-4.5M of 5.2M variants
```

**Tools:** bcftools 1.x (installed via Homebrew)

### Note on normalization

Variant normalization (`bcftools norm`) was not performed because:
1. A local GRCh37 reference FASTA was not available
2. SNPs (98%+ of variants) are unaffected by normalization
3. Only complex indels may have inconsistent representation vs dbSNP

If a GRCh37 FASTA becomes available, normalize before annotation:
```bash
bcftools norm -m -any -f hs37d5.fa genome.vcf.gz -Oz -o genome.norm.vcf.gz
tabix -p vcf genome.norm.vcf.gz
# Then annotate genome.norm.vcf.gz instead
```

## 4. Python Environment

```bash
pip3 install cyvcf2 biothings-client --break-system-packages
```

- **cyvcf2 0.31.4:** Cython/htslib VCF parser (6.9x faster than pysam)
- **biothings-client 0.4.1:** Python client for MyVariant.info API (ClinVar, gnomAD, CADD, SIFT, PolyPhen, GWAS Catalog)

## 5. MCP Server — BioMCP

Added to `~/.claude.json` mcpServers:

```json
"biomcp": {
  "type": "stdio",
  "command": "/Users/nathaYOUR_USERNAME/.local/bin/uv",
  "args": ["run", "--with", "biomcp-python", "biomcp", "run"],
  "env": {}
}
```

Provides 24 tools for variant annotation, gene info, PubMed search, and clinical trial search via MyVariant.info/MyGene.info APIs.

## 6. GRCh38 Re-Analysis Pipeline (In Progress)

### Overview

Re-analysis of the Dante Labs WGS data using modern tools to resolve known alignment artifacts
(clustered frameshifts in SON, SCARF2, AGAP3, PRKDC) and produce a clinical-grade variant report.

### Architecture

```
FASTQ (Dante Labs) -> BWA-MEM2 (GRCh38) -> DeepVariant -> VEP annotation -> Artifact triage -> Clinical report
```

Infrastructure: Imperial HPC CX3 (PBS Pro + Apptainer) for compute; local Mac for reporting.

### Pipeline: nf-core/raredisease v2.6.0

- **Aligner:** BWA-MEM2 to GRCh38 (hs38DH)
- **Variant caller:** DeepVariant v1.6+ (CNN-based, artifact-resistant)
- **SV caller:** Manta
- **Annotation:** Ensembl VEP v112+ with plugins:
  - LOFTEE (HC/LC LoF classification)
  - SpliceAI (splice disruption prediction)
  - CADD v1.7 (combined deleteriousness)
  - AlphaMissense (protein structure-based pathogenicity)
  - REVEL (ensemble missense predictor)
  - pext (proportion expressed across transcripts)
  - ClinVar
  - dbNSFP (SIFT, PolyPhen, etc.)
- **Transcript:** MANE Select prioritization

### Artifact-Aware Triage

Post-pipeline triage using slivar + vcfanno:

1. **vcfanno** adds: segmental duplications, low-complexity regions (RepeatMasker, TRF), simple repeats, gnomAD v4 frequencies
2. **slivar** flags: indel clustering (CLUSTER_ARTIFACT), allele balance skew (AB_SKEW), low quality (LOW_QUAL), biological implausibility (hom LoF in constrained genes)
3. **Tiering:**
   - Tier 1: ClinVar P/LP with 2+ review stars
   - Tier 2: LOFTEE-HC LoF in constrained genes + pext>0.1, SpliceAI>0.8, AlphaMissense pathogenic + REVEL>0.75 + CADD>25
   - Tier 3: Moderate evidence (CADD>20 + REVEL>0.5 in disease genes, moderate SpliceAI, ClinVar VUS in constrained genes)

### Kosmos Artifact Comparison

Confirms resolution of known GRCh37 alignment artifacts:
- SON (chr21): clustered frameshifts in exon 3
- SCARF2 (chr22): clustered frameshifts
- AGAP3 (chr7): clustered indels
- PRKDC (chr8): frameshift artifacts

### Read-Backed Validation

igv-reports HTML for Tier 1-2 variants:
- Mini-BAM extraction (500bp flanking)
- Bidirectional strand support
- No flanking mismatch clusters
- >=3 alt reads, MAPQ>0

### Clinical Report

Sections: QC summary, ClinVar findings, oligodontia panel (8 genes), pigmentation genetics (10 genes), high-impact non-ClinVar, LoF burden, ACMG v3.3 secondary findings (81 genes).

### HPC Configuration

```
Host: Imperial HPC CX3 (login.cx3.hpc.ic.ac.uk)
User: YOUR_USERNAME
Scheduler: PBS Pro 2024.1.1 (qsub/qstat)
Containers: Apptainer 1.4.5 (singularity CLI compatible)
Modules: Lmod 8.7.55 (tools/prod, tools/bioinf)
Key modules: Java/17, Nextflow/25.10.2, SAMtools/1.21, BCFtools/1.21
Home: /rds/general/user/YOUR_USERNAME/home/GenomeAnalysis/ (data + results)
Ephemeral: /rds/general/user/YOUR_USERNAME/ephemeral/ (Nextflow work + containers)
Config: hpc/nextflow.config (pbspro executor)
```

### Status

Phase 1 complete (2026-03-05). HPC environment verified:
- PBS job submission working (test job 1842048)
- Nextflow 25.10.2 available via module
- Apptainer 1.4.5 for container execution
- Directory structure created on RDS ephemeral
- .bashrc configured with module loading + env vars

Phase 2 in progress (2026-03-05): BAM downloaded from Dante portal, FASTQ extraction pending.

### Data Provenance (Phase 2)

BAM downloaded directly from Dante Omics portal (genome.danteomics.com) using a presigned
Cloudflare R2 URL. The portal provides time-limited download links (~1 hour validity).

- **Profile:** YOUR_PROFILE_UUID
- **File:** YOUR_SAMPLE.bam (GRCh37 alignment, ~100-150GB)
- **Download date:** 2026-03-05
- **Download method:** wget with nohup to HPC home storage

FASTQs extracted from BAM using `samtools collate | samtools fastq` (PBS job, 8 cores, 32GB RAM).
This discards the GRCh37 alignment and recovers the original paired reads for re-alignment to GRCh38.

### HPC Storage

All persistent data in home (`/rds/general/user/YOUR_USERNAME/home/GenomeAnalysis/`):
- `downloads/` — raw BAM from Dante (temporary, delete after FASTQ extraction)
- `fastq/` — pipeline-ready paired FASTQs
- `references/` — GRCh38 reference genome + annotation databases
- `results/` — nf-core/raredisease outputs

Ephemeral storage for disposable data:
- `/rds/general/user/YOUR_USERNAME/ephemeral/nf-work/` — Nextflow work directory
- `/rds/general/user/YOUR_USERNAME/ephemeral/singularity_cache/` — container images

## 7. Genome Query Skill

Created `/genome` skill at `~/.claude/skills/genome/` with:

- `SKILL.md` — Skill definition with command reference
- `genome_query.py` — Python helper (cyvcf2 + biothings-client)
- `panels/hirisplex.json` — 41-SNP HIrisPlex-S panel
- `panels/mc1r.json` — 8 MC1R red-hair variants

Trigger rules at `~/.claude/rules/genome-triggers.md`.

## 7. File Inventory

| File | Size | Description |
|------|------|-------------|
| `genome.vcf.gz` | ~1.3GB | Original VCF (backup, no rsIDs) |
| `genome.vcf.gz.tbi` | ~2MB | Tabix index for original |
| `genome.rsid.vcf.gz` | ~1.3GB | rsID-annotated VCF (primary) |
| `genome.rsid.vcf.gz.tbi` | ~2MB | Tabix index for annotated |
| `00-All.vcf.gz` | ~15GB | dbSNP b151 full reference |
| `00-All.vcf.gz.tbi` | ~3MB | dbSNP tabix index |

## 8. Software Versions

| Tool | Version | Source |
|------|---------|--------|
| bcftools | (system) | Homebrew |
| tabix/htslib | (system) | Homebrew |
| cyvcf2 | 0.31.4 | pip |
| biothings-client | 0.4.1 | pip |
| BioMCP | (latest via uv) | biomcp-python |
| Python | 3.14 | System |
| uv | 0.9.27 | ~/.local/bin |
