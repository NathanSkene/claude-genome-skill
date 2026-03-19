# Genome Query Skill Triggers

Auto-loaded rules for detecting when to invoke the /genome skill.
Copy this file to `~/.claude/rules/genome-triggers.md` to enable auto-detection.

---

## Skill Location

Main skill file: `~/.claude/skills/genome/SKILL.md`
Helper script: `~/.claude/skills/genome/genome_query.py`
Panels: `~/.claude/skills/genome/panels/`

---

## Explicit Invocation

These should invoke immediately:

| Command | Action |
|---------|--------|
| `/genome variant <rsID>` | Look up variant genotype + annotation |
| `/genome region <gene>` | Extract all variants in a gene |
| `/genome region <chr:start-end>` | Extract variants by coordinates |
| `/genome trait "<trait>"` | GWAS Catalog trait lookup |
| `/genome clinvar` | Scan for pathogenic variants |
| `/genome hirisplex` | HIrisPlex-S 41-SNP panel |
| `/genome lof` | Loss-of-function variant scan |
| `/genome impact <gene>` | Gene variants with effect scores |
| `/genome stats` | VCF statistics |
| `/genome prs` | Top/bottom polygenic risk scores |
| `/genome prs <query>` | Search PRS by trait/category/PGS ID |
| `/genome prs detail <PGS_ID>` | Detailed single PRS view |

---

## Auto-Detection Patterns (Suggest Only)

These patterns should **suggest** the genome skill, not auto-invoke.
No auto-detection triggers that make external API calls.

### Pattern 1: rsID Mentioned

When an rsID (rs followed by digits) appears in conversation:

**Detection:** `rs\d{4,}`

**Behavior:**
- Suggest: "Want me to check your genotype? `/genome variant rsXXXXXXX`"
- Do NOT auto-invoke

### Pattern 2: After /gene Disease Results

When `/gene` skill returns disease associations for a gene:

**Behavior:**
- Suggest: "Check your genome for variants in this gene? `/genome impact GENE`"
- Do NOT auto-invoke

### Pattern 3: Genotype Questions

**Trigger phrases:**
- "what's my genotype for"
- "do I carry"
- "am I a carrier"
- "check my genome"
- "look up my DNA"
- "what does my genome say"
- "my VCF"

**Behavior:**
- Invoke `/genome` with appropriate subcommand

### Pattern 4: Trait/Phenotype Questions

**Trigger phrases:**
- "hair colour" / "hair color"
- "eye colour" / "eye color"
- "skin colour" / "skin color"
- "red hair"
- "HIrisPlex"

**Behavior:**
- Suggest: "I can check the HIrisPlex-S panel for pigmentation prediction. `/genome hirisplex`"

### Pattern 5: Clinical Variant Questions

**Trigger phrases:**
- "pathogenic variants"
- "ClinVar"
- "disease risk"
- "clinical significance"
- "pharmacogenomics"

**Behavior:**
- Suggest relevant command (`/genome clinvar` or `/genome variant`)

### Pattern 6: Polygenic Risk Score Questions

**Trigger phrases:**
- "polygenic risk"
- "polygenic score"
- "PRS"
- "PGS"
- "genetic risk for"
- "my risk of"
- "risk score"
- "what's my percentile"
- "how do I compare genetically"

**Behavior:**
- Invoke `/genome prs` with appropriate subcommand
- If specific trait mentioned: `/genome prs <trait>`
- If specific PGS ID mentioned: `/genome prs detail <PGS_ID>`

### Pattern 7: Trait-Specific PRS Questions

**Trigger phrases:**
- "am I genetically predisposed to"
- "genetic risk for [trait]"
- "genetic predisposition"
- "what does my genome say about [trait]"

**Behavior:**
- Suggest: "I can check your polygenic risk score for that. `/genome prs <trait>`"
- Map common names to search terms (e.g., "heart disease" -> "coronary")

---

## VCF File Context

- **Primary VCF:** `genome.rsid.vcf.gz` (rsID-annotated)
- **Backup VCF:** `genome.vcf.gz` (original, no rsIDs)
- **Supported builds:** GRCh37/hg19, GRCh38/hg38
- **Chromosome naming:** Auto-detected (plain numbers or chr-prefixed)

---

## Privacy Reminders

When using annotation functions that call external APIs:
- `annotate_variant()` sends rsIDs to MyVariant.info
- `scan_clinvar()` sends rsIDs to MyVariant.info
- `scan_lof()` sends rsIDs to MyVariant.info
- `get_gene_impact()` sends rsIDs to MyVariant.info + Ensembl REST

Local-only operations (no external calls):
- `get_variant_by_rsid()` — queries local VCF
- `get_variant_by_pos()` — queries local VCF
- `get_region()` — queries local VCF
- `extract_panel()` — queries local VCF
