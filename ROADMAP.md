# Genome Analysis System — Roadmap

**Created:** 2026-03-19
**Last updated:** 2026-03-19
**Audience:** Normal people understanding their own genome (DTC/personal genomics)

---

## 1. What We Have Today

### nf-core/raredisease (upstream pipeline)

What it does for us:
- Alignment (BWA-MEM2)
- SNV calling (DeepVariant)
- SV calling (Manta, TIDDIT)
- CNV calling (CNVnator)
- Annotation (VEP with CADD, REVEL, SpliceAI, gnomAD, ClinVar, MaxEntScan)
- MT variant calling (Mutect2)
- Repeat expansions (ExpansionHunter, Stranger)
- ROH detection (bcftools)
- Ranking (GENMOD)
- Sample QC (Peddy)

**Status:** Pipeline completed successfully 14-Mar-2026 on Imperial CX3. 1d 13h runtime, 316.4 CPU hours. VEP annotation run standalone with 5 plugins (CADD, LOFTEE, SpliceAI, AlphaMissense, REVEL) — initial whole-genome job killed at 8h walltime, now re-running with per-chromosome parallelization (24 × 4h jobs).

### Custom features built on top

| Feature | Status | Assessment |
|---------|--------|------------|
| **Variant triage engine** (vcfanno + slivar + tiering) | Scripts ready, not yet run on raredisease output | Core value — classifies variants into actionable tiers |
| **Interactive query engine** (genome_query.py) | Working on GRCh37 VCF | Needs updating for GRCh38 output |
| **SNP panels** (HIrisPlex-S, MC1R, ACMG v3.3, pigmentation, oligodontia) | Working | Good consumer-facing features |
| **ClinVar genome-wide scan** | Working | High value for consumers |
| **LoF variant scan** | Working | Good, needs better phenotype context |
| **GWAS trait lookups** | Working | Useful but limited without PRS context |
| **IGV validation** (mini-BAM + igv-reports) | Scripts ready, not run | Important for QC but more clinical than consumer |
| **Clinical report generator** | Scripts ready, not run | Needs redesign for consumer audience |
| **HPC orchestration** (PBS Pro) | Complete, working | Infrastructure — done |
| **PRS computation** | Partially complete (121/~2900 scores) | Separate workstream |

### What hasn't been tested yet

Phases 4-6 (triage, validation, report) have scripts written but have **never been run against actual raredisease output**. This is the immediate next step before any new feature work.

---

## 2. Shortfalls in nf-core/raredisease

### What raredisease does well
- Alignment and variant calling are gold-standard (DeepVariant, GATK)
- VEP annotation with good plugin set
- Handles MT, repeat expansions, SVs
- Reproducible, containerized, well-maintained

### What raredisease does NOT do
- **No phenotype-to-variant mapping** — annotates variants but doesn't connect them to human-readable health outcomes
- **No consumer-facing interpretation** — output is a VCF, not "what does this mean for me"
- **No pharmacogenomics** — doesn't call star alleles for drug metabolism genes
- **No carrier screening** — doesn't flag recessive carrier status
- **No periodic re-analysis** — one-shot pipeline, no mechanism to check for ClinVar reclassifications over time
- **No ancestry-informed filtering** — uses global gnomAD AF, not population-specific
- **No SV interpretation** — calls SVs but doesn't interpret them (dosage sensitivity, gene overlap)
- **Singleton-only in our setup** — we run it as singleton (no trio/family analysis)

### What raredisease does poorly for our use case
- **GENMOD ranking is crude** — designed for pedigree-aware rare disease, not consumer interpretation
- **VEP output is dense and technical** — no translation layer to plain language
- **No prioritization without phenotype input** — everything is equally ranked without clinical context

---

## 3. Process Experience with nf-core/raredisease

### What went well
- Pipeline ran successfully first time after config tuning
- PBS Pro integration worked (with Apptainer)
- 1.5 day runtime for WGS is reasonable
- Resume/retry worked when individual jobs failed

### Pain points
- **Initial setup was complex** — reference staging, container caching, PBS config all needed manual work
- **Memory escalation needed tuning** — some VEP jobs needed >64GB
- **Output structure is opaque** — lots of subdirectories, not obvious which VCF to use downstream
- **GRCh38 only** — our original Dante Labs data was GRCh37, needed BAM→FASTQ→realign

### Operational notes
- Results at `/rds/general/user/nskene/home/GenomeAnalysis/results/` on CX3
- Need to transfer annotated VCF to local for downstream analysis
- SSH ControlMaster helps but CX3 doesn't accept publickey auth

---

## 4. Feature Roadmap — Ranked by Value for Consumers

### Priority 1: HPO-Driven Variant-to-Phenotype Mapping (TOP PRIORITY)

**What:** For every interesting variant found (pathogenic, likely pathogenic, LoF, high CADD, unusual), map the affected gene to HPO phenotypes and present in plain language.

**Why #1:** This is the bridge between "you have a variant in SCN1A" and "this gene is associated with febrile seizures, epileptic encephalopathy, and migraine." Without this, all other features produce opaque variant lists.

**How it works:**
1. Find variants that are unusual/deleterious/pathogenic (triage engine already does this)
2. Map affected genes → HPO phenotypes (via OMIM, Orphanet, HPO gene-disease annotations)
3. Present: gene name, variant effect, associated phenotypes in plain English, inheritance pattern, penetrance notes where available

**Data sources:**
- HPO gene-to-phenotype annotations: https://hpo.jax.org/data/annotations
- OMIM gene-disease mappings
- Orphanet gene-disease associations
- ClinGen gene-disease validity

**Difficulty:** Medium
**Dependencies:** Triage engine must work first (phases 4-6)

### Priority 2: Pharmacogenomics (PGx) Panel

**What:** Star allele calling for drug metabolism genes. Tells you how you metabolize common drugs.

**Why:** Immediately actionable. "You're a CYP2D6 poor metabolizer — codeine won't work for you, tramadol won't work, consider alternatives." Affects ~25% of prescribed drugs.

**Key genes:** CYP2D6, CYP2C19, CYP2C9, CYP3A5, DPYD, TPMT, UGT1A1, SLCO1B1, NUDT15, VKORC1

**Tools:** PyPGx or Stargazer for star allele resolution from WGS data. Map to CPIC guidelines for drug-gene pairs.

**Difficulty:** Medium-Hard (CYP2D6 is structurally complex — gene deletions, duplications, hybrid alleles)

### Priority 3: Carrier Screening Panel

**What:** Flag heterozygous pathogenic variants in recessive disease genes. "You're a carrier for cystic fibrosis."

**Why:** Reproductive planning is a primary consumer use case. ~300 genes from ACMG carrier screening recommendations.

**How:** Gene list + ClinVar P/LP filter + heterozygous genotype check. Report condition, severity, carrier frequency by ancestry.

**Difficulty:** Easy-Medium (gene list curation is the main work)

### Priority 4: Periodic Variant Re-Analysis

**What:** Re-run the same VCF against updated ClinVar annually. Flag reclassifications.

**Why:** ClinVar reclassifies ~5-8% of variants per year. A VUS today may be pathogenic next year.

**What to detect:**
- VUS promoted to P/LP
- P/LP downgraded to VUS/benign
- New ClinVar entries matching your variants
- gnomAD frequency changes with new populations

**Difficulty:** Easy (diff operation on annotation databases)

### Priority 5: Ancestry-Informed Allele Frequency Filtering

**What:** Infer genetic ancestry from WGS, use population-specific gnomAD AFs instead of global.

**Why:** A variant at 0.01% globally but 2% in your population is not rare for you. Reduces false positives.

**Difficulty:** Easy-Medium (gnomAD already has per-population AFs)

### Priority 6: SV Interpretation

**What:** Annotate structural variants against ClinGen dosage sensitivity, DGV, gnomAD-SV. Tier SVs like SNVs.

**Why:** SVs account for ~15% of pathogenic variation. Currently called but not interpreted.

**Tool:** AnnotSV

**Difficulty:** Medium-Hard

### Priority 7: MT Heteroplasmy Interpretation

**What:** Extract heteroplasmy fraction from Mutect2 output. Annotate against MITOMAP pathogenicity thresholds. Assign haplogroup via Haplogrep3.

**Why:** Many MT variants only matter above tissue-specific thresholds (>60-80%). Without this, false positives.

**Difficulty:** Easy (data already in VCF)

### Priority 8: Repeat Expansion Deep Interpretation

**What:** Map ExpansionHunter repeat counts to normal/pre-mutation/full mutation ranges for ~30 clinically significant loci.

**Why:** "You have 28 CAG repeats in HTT — this is in the normal range (below 36)" is useful. Currently just raw repeat counts.

**Difficulty:** Easy (lookup table)

### Priority 9: ACMG/AMP Auto-Classification

**What:** For variants in ACMG v3.2 secondary findings genes, auto-apply evidence codes (PS1, PM1, PM2, PP3).

**Why:** Moves from "variant list" to "classified variant list" with explicit evidence. Gets 80% of the way to expert review.

**Tools:** InterVar, AutoPVS1

**Difficulty:** Medium

### Priority 10: Family-Aware Analysis (Trio/Duo)

**What:** Support trios (proband + parents) for de novo calling, compound het detection, phasing.

**Why:** Increases diagnostic yield 15-25%. But requires family members to also be sequenced, so lower priority for single-person consumer use.

**Difficulty:** Medium-Hard

---

## 5. Immediate Next Steps

1. **~~Run VEP annotation~~** — Per-chromosome parallel VEP submitted 19-Mar-2026 (24 jobs × 4h). Scripts: `run_vep_parallel.sh` → `merge_vep.sh`
2. **Transfer complete VEP VCF from HPC to local** — after merge completes
3. **Run phases 4-6 on actual raredisease output** — triage, validation, report against the GRCh38 VCF
4. **Assess what works and what breaks** — the scripts were written before raredisease output was available
5. **Update genome_query.py for GRCh38** — currently works on GRCh37
6. **Then start Priority 1 (HPO mapping)** — this is the foundation everything else builds on

---

## 6. Architecture Notes

### Consumer vs Clinical Focus

This system is designed for **normal people understanding their own genome**, not clinical diagnostics labs. This means:
- Plain language over ACMG codes
- "What does this mean for me" over "variant classification"
- Pharmacogenomics and carrier screening over diagnostic yield
- HPO phenotype descriptions over HPO-driven filtering (we show what genes do, not filter by symptoms)

### Commercial Platform Comparison (for reference)

| Platform | What it does | Our overlap |
|----------|-------------|-------------|
| **Franklin** (Genoox) | SaaS variant classification for clinical labs. Auto-applies ACMG/AMP evidence codes. | Our ACMG auto-classification (Priority 9) |
| **Emedgene** (Illumina) | Clinical interpretation + periodic re-analysis. HPO-driven prioritization for diagnostics. | Our re-analysis (Priority 4) |
| **Fabric Genomics** | Clinical decision support for genome interpretation. Phenotype-driven. | Our HPO mapping (Priority 1) |
| **23andMe / Nebula** | Consumer genomics. Trait reports, carrier screening, PGx (limited). | Our Priorities 1-5 but deeper |

Our system aims to be deeper than 23andMe/Nebula (full WGS, not genotyping array) while being consumer-oriented rather than clinical like Franklin/Emedgene.
