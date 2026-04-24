<img align="center" src="/docs/taco-icon.png">

# TACO

**Telomere-Aware Contig Optimization**

TACO is a telomere-aware all-in-one multi-assembler comparison and refinement pipeline for genome assembly benchmarking, decision-making, and chromosome-end improvement. Developed for small eukaryotic genomes with a focus on fungal genomes, TACO runs multiple assemblers, standardizes their outputs, evaluates assembly quality, detects telomere-supported contigs, and can either **(1)** stop after generating a unified assembler comparison table for benchmarking, or **(2)** continue into telomere-aware backbone refinement for an improved chromosome-scale candidate assembly.

TACO was developed at the **Grainger Bioinformatics Center, Field Museum of Natural History**.

**What TACO is:** TACO compares multiple long-read assemblies generated from a single dataset, selects the best backbone assembly, and conservatively improves it using telomere-supported contigs from all assemblers. It produces a primary-style chromosome-level assembly with full provenance tracking and coverage QC.

**What TACO is not:** TACO does not perform Hi-C scaffolding (no YaHS/3D-DNA), does not require Hi-C data, and does not attempt full chromosome scaffolding. It is not a diploid/polyploid phasing tool. For scaffolding, use TACO's output as input to a dedicated scaffolder.

![Latest Version](https://img.shields.io/github/v/tag/yksun/TACO?sort=semver&label=Latest%20Version)
![Last Commit](https://img.shields.io/github/last-commit/yksun/TACO)
![Issues](https://img.shields.io/github/issues/yksun/TACO)
![License](https://img.shields.io/github/license/yksun/TACO)

---

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Sequencing Platform Support](#sequencing-platform-support)
- [Pipeline Steps](#pipeline-steps)
- [Telomere Detection](#telomere-detection)
- [Assembly Selection Strategy](#assembly-selection-strategy)
- [Output Structure](#output-structure)
- [Project Structure](#project-structure)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [License](#license)

---

## Overview

Genome assemblers often produce different results from the same long-read dataset. One assembler may recover longer contigs, another may preserve more complete chromosome ends, and another may provide a better balance of completeness, contiguity, and redundancy. TACO makes these comparisons systematic, interpretable, and reproducible.

TACO operates in two modes. In **assembly-only mode** (`--assembly-only`), the pipeline runs all assemblers, standardizes their outputs, evaluates quality with BUSCO, QUAST, telomere detection, and optional Merqury, then produces a unified comparison table at `assemblies/assembly_info.csv`. In **full refinement mode**, TACO continues from the comparison step into telomere-pool construction and backbone refinement, producing an improved chromosome-scale candidate assembly with preserved telomeric ends.

## Features

- Runs six long-read assemblers (HiCanu, NextDenovo, Peregrine, IPA, Flye, Hifiasm) from a single command
- Supports PacBio HiFi, Oxford Nanopore, and PacBio CLR reads via `--platform`
- Standardizes assembly outputs for direct cross-assembler comparison
- Hybrid telomere detection with de novo k-mer discovery, built-in motif families, and per-end composite scoring
- Three-tier telomere classification: strict T2T, single-end strong, and telomere-supported
- Benchmarks assemblies with BUSCO, QUAST, telomere metrics, and optional Merqury
- Biologically informed automatic backbone selection (smart scoring)
- Assembly-only mode for convenient benchmarking without refinement
- Telomere-aware backbone refinement with redundancy reduction and telomeric end rescue
- Machine-readable benchmark logs for reproducible reporting

## Installation

See [INSTALLATION.md](INSTALLATION.md) for detailed instructions.

### Quick Install

```bash
git clone https://github.com/yksun/TACO.git
cd TACO
conda env create -f taco-env.yml
conda activate taco
pip install -e .

# Verify
taco --help
```

After installation, the `taco` command is available system-wide within the conda environment.

### Requirements

TACO requires a Unix-like system (Linux or macOS) with Python >= 3.8 and Conda. Most dependencies are installed via `taco-env.yml`. TACO's Python modules use only the standard library with no additional pip packages. Canu, Peregrine, and IPA require manual installation (see [INSTALLATION.md](INSTALLATION.md)). If any assembler is absent or fails, TACO skips that step and continues with the others.

## Quick Start

**Full refinement run** (comparison + telomere-aware refinement):

```bash
mkdir -p my_project && cd my_project

taco -g 12m -t 16 \
  --fastq /path/to/reads.fastq \
  --taxon fungal
```

**Assembly-only comparison** (benchmarking only):

```bash
taco -g 12m -t 16 \
  --fastq /path/to/reads.fastq \
  --taxon fungal \
  --assembly-only
```

**With Merqury QV scoring (optional, auto-detected if installed):**

Merqury is automatically enabled when `merqury.sh` is on PATH and a `.meryl` database is found in the working directory. You can also specify the database explicitly:

```bash
taco -g 12m -t 16 \
  --fastq /path/to/reads.fastq \
  --taxon fungal \
  --merqury-db reads.meryl
```

**With explicit motif override** (only if biologically known):

```bash
taco -g 12m -t 16 \
  --fastq /path/to/reads.fastq \
  -m TTAGGG
```

## Usage

### Parameters

| Parameter | Description |
|---|---|
| `-g`, `--genomesize` | Estimated haploid genome size (e.g., `12m`, `40m`, `2g`) |
| `-t`, `--threads` | Number of CPU threads |
| `--fastq` | Input FASTQ file (use absolute path) |
| `--taxon` | Taxonomy preset for telomere detection: `vertebrate`, `animal`, `plant`, `insect`, `fungal`, or `other` (default). Sets motif-family priors and detection behavior automatically. |
| `-m`, `--motif` | Telomere motif override (optional). Only use when the exact motif is biologically known for the species. When omitted, taxon-aware hybrid detection is used instead. |
| `--platform` | Sequencing platform: `pacbio-hifi` (default), `nanopore`, or `pacbio`. Also determines default polishing tool. |
| `-s`, `--steps` | Run selected steps only (e.g., `1,3-5`) |
| `--reference`, `-ref` | Reference FASTA for comparison. Included as the "reference" assembler in all comparison tables. |
| `--busco` | Run BUSCO (optionally specify lineage dataset) |
| `--choose` | Manually choose the backbone assembler |
| `--assembly-only` | Stop after assembler comparison |
| `--auto-mode` | Backbone selection mode: `smart` (default) or `n50` |
| `--merqury` | Force-enable Merqury (auto-detected if `merqury.sh` + `.meryl` db found) |
| `--merqury-db` | Enable Merqury with a specific `.meryl` database path |
| `--no-merqury` | Disable Merqury even if installed and auto-detected |
| `--no-purge-dups` | Skip purge_dups after refinement |
| `--no-polish` | Skip automatic polishing after refinement |
| `--allow-t2t-replace` | Allow rescue donors to replace immutable Tier 1 (T2T) contigs. Disabled by default for safety |

### Assembly-Only Mode

Use `--assembly-only` when the goal is assembler benchmarking and comparison without refinement. TACO runs all assemblers, standardizes outputs, runs BUSCO, telomere detection, QUAST, and optional Merqury, then writes the combined comparison table to `assemblies/assembly_info.csv` and a summary to `final_results/assembly_only_result.csv`.

## Sequencing Platform Support

TACO supports three sequencing platforms. Each assembler receives platform-appropriate flags automatically. The platform also determines the default polishing strategy: HiFi assemblies are polished with NextPolish2 by default (k-mer-based, safe for high-accuracy reads; requires `yak` for k-mer database construction), Nanopore assemblies use Medaka (neural-network polisher; falls back to Racon), and CLR assemblies use Racon.

| Platform | `--platform` | Assemblers Used | Notes |
|---|---|---|---|
| PacBio HiFi | `pacbio-hifi` (default) | canu, nextDenovo, peregrine, IPA, flye, hifiasm | All 6 assemblers |
| Oxford Nanopore | `nanopore` | canu, nextDenovo, flye | Peregrine, IPA, hifiasm skipped |
| PacBio CLR | `pacbio` | canu, nextDenovo, peregrine, flye | IPA, hifiasm skipped |

Incompatible assemblers are automatically skipped with a warning. IPA and hifiasm only support PacBio HiFi reads. Peregrine does not support Nanopore reads.

### Platform-Specific Assembler Flags

Each assembler receives the appropriate read-type flag automatically:

| Assembler | HiFi | Nanopore | CLR |
|---|---|---|---|
| Canu | `-pacbio-hifi` | `-nanopore` | `-pacbio` |
| Flye | `--pacbio-hifi` | `--nano-hq` (Q20+) | `--pacbio-raw` |
| Hifiasm | default mode | skipped | skipped |
| IPA | default mode | skipped | skipped |
| Peregrine | default mode | skipped | default mode |
| NextDenovo | via config file | via config file | via config file |

For older ONT data that is not Q20+ basecalled, set `FLYE_ONT_FLAG=--nano-raw` in the environment.

### Platform-Specific Polishing Strategy

| Platform | Polishing Tool | Notes |
|---|---|---|
| HiFi | NextPolish2 (yak k-mer based) | K-mer polishing corrects residual errors safely; requires `nextpolish2` + `yak` |
| Nanopore | Medaka (Racon fallback) | Neural-network polisher; set `MEDAKA_MODEL` for non-default chemistry |
| CLR | Racon | Standard error-correction for CLR reads |

### Combined Platform × Taxon Strategy Overview

| Component | Fungal | Plant | Vertebrate/Animal | Insect/Other |
|---|---|---|---|---|
| **Default BUSCO lineage** | ascomycota_odb10 | embryophyta_odb10 | vertebrata_odb10 / metazoa_odb10 | insecta_odb10 / requires `--busco` |
| **Merqury** | auto if HiFi | auto if HiFi | auto if HiFi | auto if HiFi |
| **Telomere motifs** | TTAGGG + TG1-3 + Candida | TTTAGGG | TTAGGG | TTAGG (insect) / all (other) |
| **Score window** | 300 bp | 1000 bp | 1000 bp | 500 bp |
| **Backbone scoring** | S×1000 - D×600 + T2T×350 | S×1000 - D×300 + T2T×200 | S×1000 - D×500 + T2T×200 | S×1000 - D×500 + T2T×300 |
| **BUSCO trial C-drop** | 2% (strict) | 4% (relaxed) | 3% (moderate) | 2% (default) |
| **purge_dups mode** | two-round (haploid-aggressive) | single-round + polyploid warning | two-round | single-round |
| **Polishing (HiFi)** | NextPolish2 (yak k-mer based) | NextPolish2 (yak k-mer based) | NextPolish2 (yak k-mer based) | NextPolish2 (yak k-mer based) |
| **Polishing (ONT)** | Medaka → Racon | Medaka → Racon | Medaka → Racon | Medaka → Racon |
| **Polishing (CLR)** | Racon | Racon | Racon | Racon |

## Pipeline Steps

| Step | Description |
|---|---|
| 1 | HiCanu assembly |
| 2 | NextDenovo assembly |
| 3 | Peregrine assembly |
| 4 | IPA assembly |
| 5 | Flye assembly |
| 6 | Hifiasm assembly |
| 7 | Copy and normalize all assemblies |
| 8 | BUSCO on all assemblies |
| 9 | Telomere contig detection and telomere metrics |
| 10 | Build optimized telomere pool across assemblies |
| 11 | QUAST for assembler comparison |
| 12 | Backbone selection and telomere-aware refinement |
| 13 | BUSCO on final assembly |
| 14 | Telomere analysis of final assembly |
| 15 | QUAST on final assembly |
| 16 | Final comparison report |
| 17 | Cleanup into structured output folders |
| 18 | Assembly-only comparison summary |

With `--assembly-only`, TACO runs Steps 1-9, 11, and 18 (skips Step 10 which builds the refinement telomere pool) and stops before backbone refinement.

## Telomere Detection

TACO v1.3.0 uses a taxon-aware hybrid telomere detection system that combines built-in motif families with de novo k-mer discovery.

### Taxon-Aware Presets

Use `--taxon` to select the appropriate telomere motif priors for your organism. This is the recommended approach instead of forcing `--motif` directly.

| `--taxon` | Primary Motifs | Notes |
|---|---|---|
| `vertebrate` | TTAGGG | Highly conserved; exact motif matching is most reliable here |
| `animal` | TTAGGG | Strong prior for vertebrates, less certain for distant metazoans |
| `plant` | TTTAGGG | Common plant repeat; some lineages vary |
| `insect` | TTAGG | Common insect repeat; not universal across all insect orders |
| `fungal` | TTAGGG, TG1-3, Candida | Diverse fungal telomeres — all built-in families used |
| `other` (default) | All families | Unknown taxon — relies on de novo discovery plus all priors |

The `--motif` flag is an optional override. Do not force `--motif` unless the telomere repeat is biologically confirmed for your species or lineage. For fungi and unknown taxa especially, forcing a motif may miss true telomeres.

### Built-in Motif Families

TACO ships with five motif families: the canonical vertebrate/filamentous fungal TTAGGG repeat, the budding yeast TG1-3/C1-3A degenerate repeat, the Candida 23-bp repeat (ACGGATGTCTAACTTCTTGGTGT), the plant TTTAGGG repeat, and the insect TTAGG repeat.

### Scoring System

Each contig end is scored using a composite of four metrics: telomere density (weight 0.40), longest consecutive run (weight 0.30), distance of repeats from the contig terminus (weight 0.20), and covered base pairs (weight 0.10). The composite score ranges from 0 to 1.

### Classification Tiers

Contigs are classified into three tiers based on their end scores: **strict T2T** contigs have strong telomere signal at both ends (score >= 0.25 at each end), **single-end strong** contigs have strong signal at one end only, and **telomere-supported** contigs have at least weak signal (score >= 0.08) at one end.

## Assembly Selection Strategy

When `--choose` is not provided, TACO automatically selects the backbone assembly for refinement. The scoring formula adapts its weights based on `--taxon` to match the biological characteristics of each organism group.

### Smart Scoring (default)

TACO ranks assemblies using a composite score. Merqury QV and completeness are included when available (optional — auto-detected if `merqury.sh` is installed and a `.meryl` database is found; otherwise these terms contribute 0):

```
score = BUSCO_S × w_busco_s + T2T × w_t2t + single_tel × w_single
      + MerquryComp × 200 + MerquryQV × 20      ← optional, 0 if Merqury not available
      - contigs × w_contigs + log10(N50) × w_n50
      - BUSCO_D × w_busco_d
```

BUSCO single-copy completeness (S%) is used instead of total completeness (C%) to avoid rewarding highly duplicated assemblies. BUSCO duplication (D%) is explicitly penalised. When Merqury is available, its k-mer-based QV and completeness provide an independent quality signal that helps distinguish assemblies with similar BUSCO scores. The weights are tuned per taxon as described below.

### Taxon-Specific Scoring Strategies

**Fungal** (`--taxon fungal`): Fungal genomes are typically small (10–60 Mb) with well-defined chromosomes. TACO uses strict BUSCO duplicate penalty (`w_busco_d = 600`) because duplicated assemblies are almost always artefactual in haploid fungi. T2T contigs are weighted heavily (`w_t2t = 350`) since telomere rescue is highly effective for small genomes where individual T2T chromosomes can be resolved. The contig-count penalty remains moderate (`w_contigs = 30`) because most fungal genomes have few chromosomes.

**Plant** (`--taxon plant`): Plant genomes vary enormously in size and ploidy. TACO relaxes the BUSCO duplicate penalty (`w_busco_d = 300`) because polyploidy naturally inflates D% even in correct assemblies. The contig-count penalty is increased (`w_contigs = 50`) to discourage fragmented assemblies in these often large genomes. T2T weight is reduced (`w_t2t = 200`) because long repetitive arrays near telomeres can produce false-positive signals, and interstitial telomeric repeats (ITRs) are common in plants.

**Vertebrate / Animal** (`--taxon vertebrate` or `--taxon animal`): Vertebrate genomes are large (1–3+ Gb) and repeat-rich. TACO increases the N50 weight (`w_n50 = 200`) to favour contiguous assemblies and moderates the contig-count penalty (`w_contigs = 40`). T2T weight is reduced (`w_t2t = 200`) because interstitial telomeric repeats are frequent in vertebrates and can inflate telomere counts. BUSCO duplicate penalty stays at the default (`w_busco_d = 500`).

**Insect / Other** (`--taxon insect` or `--taxon other`): These taxa use the balanced default weights: `w_busco_s = 1000`, `w_t2t = 300`, `w_single = 150`, `w_contigs = 30`, `w_n50 = 150`, `w_busco_d = 500`. This is appropriate when the biological characteristics of the target organism are not well characterised.

| Weight | Fungal | Plant | Vertebrate/Animal | Insect/Other |
|---|---|---|---|---|
| `w_busco_s` | 1000 | 1000 | 1000 | 1000 |
| `w_t2t` | 350 | 200 | 200 | 300 |
| `w_single` | 150 | 150 | 150 | 150 |
| `w_contigs` | 30 | 50 | 40 | 30 |
| `w_n50` | 150 | 150 | 200 | 150 |
| `w_busco_d` | 600 | 300 | 500 | 500 |

### Taxon-Specific Telomere Detection Windows

The default telomere score window also varies by taxon to match typical telomere array lengths: fungi use 300 bp (fungal telomere arrays are short, often 50–300 bp), plants and vertebrates use 1000 bp (longer repeat arrays), and other taxa use 500 bp (balanced default).

### Taxon-Specific BUSCO Trial Thresholds (Step 12F)

When validating rescue candidates via BUSCO trial, the maximum acceptable C% drop, M% rise, and D% rise depend on taxon:

- **Fungi**: strict thresholds (2% C-drop, 0.3% M-rise, 2% D-rise) — haploid genomes with stable BUSCO profiles.
- **Plant**: relaxed thresholds (4% C-drop, 1.0% M-rise, 6% D-rise) — polyploidy causes natural BUSCO variability.
- **Vertebrate**: moderate thresholds (3% C-drop, 0.5% M-rise, 4% D-rise).
- **Other**: balanced defaults (2.5% C-drop, 0.5% M-rise, 3% D-rise).

A D-rise (duplicated BUSCO increase) check catches cases where a rescue introduces redundant copies of single-copy orthologs — a sign of retained haplotigs or mis-joined contigs. All thresholds can be overridden via `STEP12_MAX_BUSCO_C_DROP`, `STEP12_MAX_BUSCO_M_RISE`, and `STEP12_MAX_BUSCO_D_RISE` environment variables.

Additional environment variables for fine-tuning Step 12: `PROTECT_COV` / `PROTECT_ID` (strict dedup thresholds, default 0.95/0.95), `DEDUP_MAX_BUSCO_C_DROP` (maximum tolerated BUSCO C drop after dedup, default 3.0%), `CHIMERA_MIN_CROSS_COV` (minimum cross-assembly coverage for chimera mapping check, default 0.60), `AGGR_NONTELO_COV` / `AGGR_NONTELO_ID` (taxon-aware non-telo dedup), `SELFDEDUP_COV` / `SELFDEDUP_ID` (self-dedup thresholds).

### Taxon-Specific Rescue Limits (Step 12F)

The maximum number of accepted rescue candidates per run is taxon-aware: fungi allow up to 20 rescues (many small chromosomes), vertebrates 10, plants 8 (conservative due to polyploidy risk), and other taxa 15. This prevents runaway replacement in complex genomes.

### D-Aware Duplicate Filter (Step 12F2)

Before adding "novel" pool T2T contigs, TACO checks each candidate against the current backbone with minimap2 and applies a five-tier decision:

1. **Overlaps Tier 1 (T2T) backbone** → reject (pure duplicate, backbone already has T2T).
2. **Overlaps Tier 2 (non-T2T) backbone at ≥80% target coverage** → upgrade (replace backbone with T2T — better telomere evidence, comparable size).
3. **Overlaps Tier 2 at 50–80% target coverage** → read-coverage diagnostic. TACO maps HiFi reads to the backbone contig and compares median coverage in the T2T-covered region vs the uncovered region. If the uncovered region has very low coverage (< 30% of covered), the backbone is chimeric and the T2T is the real chromosome — TACO replaces the backbone. If coverage is normal, the backbone is real — TACO rejects the novel contig (adding it would increase BUSCO D). Configurable via `CHIMERIC_COV_RATIO` (default 0.30).
4. **Overlaps Tier 2 at < 50%** → reject (insufficient overlap for any useful decision).
5. **No significant overlap** → add as genuinely novel chromosomal region (with optional BUSCO D check: `NOVEL_MAX_D_RISE`).

Thresholds: `NOVEL_DUP_COV` (default 0.80), `NOVEL_DUP_ID` (default 0.90), `NOVEL_UPGRADE_TCOV` (default 0.80).

### Taxon-Specific purge_dups Behaviour (Step 12H)

purge_dups strategy is taxon-aware. Fungi and other haploid genomes use two-round purging (`-2` flag) for more aggressive duplicate detection — in haploid genomes, duplicated contigs receive similar read coverage to primary contigs, making them harder to detect with single-round purging. Vertebrate and animal genomes also use two-round purging for thorough haplotig removal. Plant genomes use conservative single-round purging to avoid collapsing homeologous sequences in polyploid species — use `--no-purge-dups` if this is still too aggressive. Coverage cutoffs from `calcuts` are logged for debugging; override with the `PURGE_DUPS_CALCUTS` environment variable if automatic thresholds are incorrect for your dataset.

### N50-only Mode

`--auto-mode n50` selects the assembly with the highest N50. This reproduces legacy behavior but may favor contiguous assemblies that lack completeness.

## Step 12 — T2T-First Telomere-Aware Backbone Refinement

Step 12 adopts a T2T-first assembly philosophy with a **two-tier confidence model**:

- **Tier 1 (Immutable)**: T2T contigs — contigs with verified telomere signal at both ends. These are treated as protected chromosomal anchors and are never replaced during rescue, unless `--allow-t2t-replace` is explicitly set. This protects the highest-confidence contigs from accidental degradation.
- **Tier 2 (Editable)**: Backbone contigs — gap-fill contigs that cover chromosomal regions not represented by T2T contigs. These may be replaced by telomere-bearing rescue donors if the replacement passes BUSCO trial validation.

Backbone contigs are preserved by default to maintain BUSCO completeness — purge_dups at Step 12H handles haplotig removal using read-coverage evidence. Novel T2T additions undergo a D-aware duplicate filter that either upgrades Tier 2 backbone contigs (replacing non-T2T with T2T) or rejects pure duplicates of Tier 1 contigs. Rescue donors must carry verified telomere signal.

### Step 12 Sub-step Flow

1. **12A** — Merqury QV scoring (optional; auto-detected if installed, or enabled with `--merqury`/`--merqury-db`).
2. **12B** — auto-select backbone assembler (smart scoring with taxon-aware weights).
3. **12C** — prepare cleaned backbone + chimera safety using two strategies: (a) **size gate** — contigs > 1.5× the largest individual assembler contig are flagged; (b) **cross-assembly mapping** — each protected contig is aligned against all other assembler outputs; contigs not well-covered (≥60%) by any single assembler's contig are flagged as potential chimeras. Configurable via `CHIMERA_MIN_CROSS_COV`.
4. **12D** — backbone-first classification and pool T2T analysis:
   - **12D1** backbone telomere classification: classify backbone contigs as Tier 1 (T2T, immutable) or Tier 2 (non-T2T, upgradeable).
   - **12D2** pool T2T analysis: align pool T2T contigs against backbone. Pool T2T redundant to Tier 1 backbone are discarded. Pool T2T that cover a Tier 2 backbone contig (≥80% target coverage, ≥85% identity) become upgrade donors. Others are candidate novel additions.
   - **12D3** backbone self-dedup: disabled by default (preserves BUSCO completeness). purge_dups at 12H handles haplotig removal. Re-enable with `SELFDEDUP_ENABLE=1`.
5. **12E** — telomere upgrade: pool T2T donors replace Tier 2 backbone contigs. Tier 1 (T2T) contigs are immutable — candidates targeting them are rejected unless `--allow-t2t-replace` is set. Each replacement is assigned a class: `upgrade_tier2_to_t2t`, `replace_single_with_better`, etc.
6. **12F** — BUSCO trial validation: for each candidate, build a trial assembly and run BUSCO. Rejection thresholds are taxon-aware (fungi: 2% C-drop / 2% D-rise, plant: 4% / 6%, vertebrate: 3% / 4%). D-aware novel filter (12F2): novel T2T contigs that overlap Tier 2 backbone REPLACE the backbone (upgrade); those overlapping Tier 1 are rejected as duplicates; those with no overlap are added as genuinely novel. Optional BUSCO D check rejects additions that increase duplication excessively.
7. **12G** — final combine: backbone (with upgrades) + novel additions. Post-upgrade dedup protects all backbone contigs; only novel additions can be removed if redundant.
8. **12H** — purge_dups: taxon-aware haplotig/duplicate purging (skip with `--no-purge-dups`).
9. **12I** — automatic polishing: NextPolish2 for HiFi (k-mer-based via yak; skip with `--no-polish`), Medaka for ONT (Racon fallback), Racon for CLR.
10. **12J** — telomere-aware genome-size pruning: only non-telomeric contigs are removed when assembly exceeds the size budget. Telomere-bearing contigs are never pruned.
11. **12K** — final assembly coverage QC: maps HiFi reads to the final assembly, computes sliding-window coverage (default 5 kb), and flags zero-coverage gaps, very-low-coverage regions, and sudden coverage drops. Outputs `coverage_summary.tsv`, `weak_regions.tsv`, and `weak_regions.gff3` (loadable in IGV).
12. **12L** — "do no harm" safety comparison: compares final assembly vs original backbone for size, telomere count, and genome size deviation. If refinement degraded quality, both assemblies are preserved with a `refinement_warning.txt` explaining the issues.

### BUSCO Trial Validation

TACO validates each telomere rescue candidate by building a trial assembly where one backbone contig is replaced by one donor contig, then running BUSCO with the same lineage selected by the user. Rejection is triggered by three independent BUSCO metrics: C% drop (completeness loss), M% rise (missing gene increase), and D% rise (duplicated BUSCO increase, catching retained haplotigs). Rejection thresholds are taxon-aware: fungi use strict thresholds (2% C-drop, 2% D-rise), plants use relaxed thresholds (4% C-drop, 6% D-rise, accounting for polyploidy), and vertebrates use moderate thresholds (3% C-drop, 4% D-rise). An additional safety check rejects `replace_single_with_better` candidates if telomere evidence weakens at either end after replacement (suspicious size drops >30% also trigger rejection). This greedy, sequential approach ensures that each accepted rescue improves or maintains assembly quality. The trial summary TSV includes `replacement_class` and `D` (duplicated %) columns for full traceability.

### Post-Refinement Stack

After the rescued/combined assembly is produced, TACO runs purge_dups by default to remove leftover haplotigs, overlapping fragments, and residual duplicates. purge_dups strategy is taxon-aware: fungi and other haploid genomes use two-round purging (`-2`) for aggressive duplicate detection — in haploid genomes, duplicated contigs receive similar read coverage to primary contigs, making them harder to detect. Vertebrate and animal genomes also use two-round purging. Plants use conservative single-round purging to preserve homeologs; a warning is emitted for polyploid risk. Coverage cutoffs from `calcuts` are logged; override with `PURGE_DUPS_CALCUTS` if automatic thresholds are incorrect. This is followed by automatic polishing selected from `--platform`: HiFi assemblies use NextPolish2 (maps HiFi reads to assembly with minimap2, sorts BAM with samtools, builds yak k-mer databases, then runs k-mer-based correction), Nanopore uses Medaka (falls back to Racon), and CLR uses Racon. Both steps can be skipped with `--no-purge-dups` and `--no-polish`.

### Diploid and Polyploid Note

TACO is designed for producing a best primary-style chromosome-level assembly, not a fully phased diploid or polyploid reconstruction. For strongly diploid or polyploid genomes, telomere-bearing contigs from different haplotypes may appear as rescue donors, and purge_dups may collapse alternative haplotigs. This is acceptable when the goal is a cleaned primary reference assembly.

### Provenance GFF3

TACO writes a GFF3 annotation file (`final.merged.provenance.gff3`) alongside the final assembly. Each contig gets one GFF3 record (type=contig) spanning its full length, with attributes documenting its full provenance chain: `source_assembler` (which assembler produced it), `assembler_contig` (the original contig name from that assembler before Step 10 pool renaming), `source_type` (assembler or quickmerge), `role` (backbone, upgrade_donor, or novel_t2t), `replacement_class` (for upgrade donors), `replaced_contig` (which backbone contig was replaced), and `description` (a human-readable summary like "Entire replacement: peregrine contig 'contig_5' replaced by canu contig 'tig00000015' (class: upgrade_tier2_to_t2t)").

For quickmerge-derived contigs, the GFF3 includes additional contig-level attributes (`qm_assembler1`, `qm_assembler2`) identifying the two source assemblers, plus child records (type=region) with `Parent` linking to the contig. Each region record spans the genomic coordinates contributed by a specific assembler, with `source_assembler` and `assembler_contig` showing the original source. For example, a quickmerge contig produced from canu × flye will have region records like "Region 1-500000 from canu contig 'tig00001'" and "Region 400000-900000 from flye contig 'contig_3'", enabling users to trace every base pair back to its assembler of origin.

A companion file `pool_contig_provenance.tsv` maps every pool contig back to its source assembler and original contig name, with extended columns for quickmerge contigs: `qm_assembler1`, `qm_assembler2`, and `qm_regions` (semicolon-delimited `start-end:assembler:contig` entries).

### Coverage QC GFF3 and Reports

TACO maps HiFi reads back to the final assembly (Step 12K) and scans for assembly errors using a sliding-window coverage analysis (default 5 kb window). Three output files are produced:

`coverage_summary.tsv` — per-contig coverage statistics: median, mean, min, max, zero-coverage bases, and low-coverage bases. Use this to identify contigs with overall poor read support.

`weak_regions.tsv` — every flagged window with precise coordinates, contig total length, source assembler, source type, window median coverage, global median coverage, coverage ratio, and a flag classifying the issue: `ZERO_GAP` (>50% of window has zero coverage — likely a misjoin), `VERY_LOW` (window median < 15% of global — possibly chimeric), `LOW` (< 30% of global — suspicious drop), `MIXED_LOW` (>30% of bases at zero or below threshold — intermittent problems).

`weak_regions.gff3` — the same weak regions as a GFF3 annotation file that can be loaded directly in genome browsers (IGV, JBrowse, etc.) alongside the final assembly FASTA. Each record has `type=coverage_warning` with the window median coverage as the score column (for color-coding by severity). Attributes include `flag`, `window_median`, `global_median`, `ratio`, `source_assembler`, `source_type`, and a human-readable `description`. Load `final.merged.fasta` as the genome and `weak_regions.gff3` as an annotation track to visually inspect problem regions.

Example workflow for inspecting weak spots:
```
# In IGV or JBrowse:
# 1. Load final_results/final.merged.fasta as genome
# 2. Load final_results/final.merged.provenance.gff3 as track (provenance)
# 3. Load final_results/weak_regions.gff3 as track (coverage warnings)
# 4. Navigate to flagged regions to inspect assembly quality
```

## Output Structure

```
project_directory/
├── assemblies/
│   ├── assembly_info.csv                # Unified comparison table
│   ├── canu.fasta                       # Normalized assembly outputs
│   ├── nextdenovo.fasta
│   ├── ...
│   ├── single_tel.replaced.debug.tsv    # All rescue alignment hits
│   ├── single_tel.candidates.tsv        # Plausible rescue candidates
│   ├── rescue_rejection_summary.txt     # Rejection reasons
│   ├── quickmerge_validation.tsv       # Quickmerge structural validation decisions
│   ├── telomere_pool_decisions.tsv     # Per-contig pool classification decisions
│   ├── rescue_trial_summary.tsv         # BUSCO trial results (with replacement_class, D%)
│   ├── final_merge.raw.fasta            # Pre-purge combined assembly
│   └── *.busco/                         # BUSCO results per assembly
├── final_results/
│   ├── final_result.csv                 # Final comparison report
│   ├── final_assembly.fasta             # Refined assembly (full mode)
│   ├── final.merged.provenance.gff3     # GFF3 provenance: full assembler tracing per contig
│   ├── pool_contig_provenance.tsv       # Pool contig → assembler + original name mapping
│   ├── coverage_summary.tsv             # Per-contig coverage stats (median, mean, zero/low bp)
│   ├── weak_regions.tsv                 # Flagged weak windows with coords, source assembler, flag
│   ├── weak_regions.gff3                # GFF3 coverage warnings: load in IGV to see weak spots
│   ├── {backbone}.backbone.original.fasta  # Original backbone (for do-no-harm comparison)
│   ├── refinement_warning.txt           # Quality warnings if refinement degraded backbone
│   └── assembly_only_result.csv         # Comparison summary (assembly-only)
├── telomere_pool/                       # Telomere pool intermediates
├── quast_results/                       # QUAST output
├── logs/                                # Per-step log files
├── benchmark_logs/                      # Machine-readable benchmark data
│   ├── run_metadata.tsv
│   ├── step_benchmark.tsv
│   └── run_summary.txt
└── version.txt                          # Software versions
```

## Project Structure

```
TACO/
├── setup.py                # pip install entry point
├── run_taco                # Shell wrapper (no install needed)
├── taco/                   # Python package
│   ├── __init__.py         # Package metadata (v1.3.0)
│   ├── __main__.py         # CLI entry point: taco [options]
│   ├── cli.py              # Argument parsing
│   ├── pipeline.py         # Pipeline runner, logging, benchmarking
│   ├── steps.py            # All 18 step implementations
│   ├── utils.py            # Shared utilities and FASTA I/O
│   ├── telomere_detect.py  # Hybrid telomere detection engine
│   ├── telomere_pool.py    # Telomere pool classification
│   ├── clustering.py       # Minimap2-based contig clustering
│   ├── backbone.py         # Backbone selection and scoring
│   └── reporting.py        # Final report generation
├── docs/                   # Documentation and images
├── taco-env.yml            # Conda environment
├── INSTALLATION.md
├── README.md
├── LICENSE
└── .gitignore
```

## Troubleshooting

**Canu reports `master +XX changes` or Step 1 fails with a Java error:** The conda environment now includes `openjdk>=11` to fix the Java runtime. If you still see this error, the bioconda canu may be a dev build — download a stable binary from https://github.com/marbl/canu/releases and place it on PATH. TACO detects dev builds and warns you; if canu fails, the pipeline continues with the remaining assemblers.

**IPA or Peregrine skipped:** These assemblers only support certain platforms. IPA requires PacBio HiFi; Peregrine does not support Nanopore. Use `--platform` to match your data type.

**Telomere motif appears incorrect:** Do not force `--motif` unless the telomere repeat is biologically known for your species. Use `--taxon` to select the appropriate preset instead. TACO's built-in motif families cover canonical TTAGGG, budding yeast TG1-3, Candida, plant TTTAGGG, and insect TTAGG repeats.

**purge_dups or polishing not running:** These tools must be installed in the conda environment. Use `conda install -c bioconda purge_dups nextpolish2 yak racon medaka` or skip with `--no-purge-dups` / `--no-polish`. For HiFi polishing, NextPolish2 and yak are both required. For Nanopore polishing, Medaka is preferred; if unavailable, Racon is used as fallback.

**`TACO.sh: command not found`:** Add the TACO directory to your PATH or run with the full path.

**Missing Python modules:** TACO uses only the Python standard library. If you see import errors, ensure Python >= 3.8 is installed and the `taco/` directory is alongside `TACO.sh`.

**Merqury not working:** Merqury is optional. Install with `conda install -c bioconda merqury meryl` or use `--no-merqury` to skip.

## Citation

If you use TACO in a publication, please cite the software and archive the exact release used for reproducibility (e.g., via Zenodo).

TACO was developed at the Grainger Bioinformatics Center, Field Museum of Natural History, Chicago, Illinois, USA.

## Changelog

See [docs/CHANGELOG.md](docs/CHANGELOG.md) for the full version history, including detailed notes on the v0.5.6 → v1.0.0 Bash-to-Python conversion.

## License

TACO is released under the [MIT License](LICENSE).
