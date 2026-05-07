<img align="center" src="/docs/taco-icon.png">

# TACO

**Telomere-Aware Contig Optimization for long-read genome assembly comparison and refinement**

TACO is a reproducible multi-assembler pipeline for long-read genome assembly benchmarking, backbone selection, and conservative chromosome-end improvement. It was developed for small to moderate eukaryotic genomes, with a particular focus on fungal genomes, but also includes taxon-aware settings for plants, vertebrates, insects, and other organisms.

From one input read set, TACO runs compatible assemblers, standardizes their outputs, evaluates assembly quality with BUSCO, QUAST, telomere metrics, and optional Merqury, then follows one of two workflows: **assembly-only benchmarking** or **full telomere-aware refinement**.

TACO was developed at the **Grainger Bioinformatics Center, Field Museum of Natural History**.

**What TACO is:** TACO compares multiple long-read assemblies generated from a single dataset, selects the best-supported backbone assembly, and conservatively improves it using telomere-supported contigs from all assemblers. It produces a primary-style chromosome-level candidate assembly with provenance tracking, final QC, and coverage diagnostics.

**What TACO is not:** TACO does not perform Hi-C scaffolding, does not require Hi-C data, and does not attempt full chromosome scaffolding. It is not a diploid/polyploid phasing tool. For scaffolding, use TACO output as input to a dedicated scaffolder such as YaHS or 3D-DNA.

![Latest Version](https://img.shields.io/github/v/tag/yksun/TACO?sort=semver&label=Latest%20Version)
![Last Commit](https://img.shields.io/github/last-commit/yksun/TACO)
![Issues](https://img.shields.io/github/issues/yksun/TACO)
![License](https://img.shields.io/github/license/yksun/TACO)

---

## Table of Contents

- [Overview](#overview)
- [Workflow At A Glance](#workflow-at-a-glance)
- [Key Features](#key-features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Command-Line Reference](#command-line-reference)
- [Sequencing Platform Support](#sequencing-platform-support)
- [Pipeline Steps](#pipeline-steps)
- [Compare Mode (`--compare`)](#compare-mode---compare)
- [Telomere Detection](#telomere-detection)
- [Assembly Selection Strategy](#assembly-selection-strategy)
- [Outputs And Reports](#outputs-and-reports)
- [Project Structure](#project-structure)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [License](#license)

---

## Overview

Genome assemblers often produce different results from the same long-read dataset. One assembler may recover longer contigs, another may preserve more complete chromosome ends, and another may provide a better balance of completeness, contiguity, and redundancy. TACO makes these comparisons systematic, interpretable, and reproducible.

TACO operates in two modes. In **assembly-only mode** (`--assembly-only`), the pipeline runs all compatible assemblers, standardizes their outputs, evaluates quality with BUSCO, QUAST, telomere detection, and optional Merqury, then produces a unified comparison table. In **full refinement mode** (default), TACO continues into telomere-pool construction, backbone refinement, final QC, and report cleanup.

## Workflow At A Glance

| Goal | Command mode | Steps run | Main result |
|---|---|---|---|
| Compare assemblers without changing any assembly | `--assembly-only` | 0-10, then 14B | `assemblies/assembly_info.csv` and `final_results/assembly_only_result.csv` |
| Produce a refined candidate assembly | default full mode | 0-14 | `final_results/final_assembly.fasta` and `final_results/final_result.csv` |
| Resume final QC and reporting after refinement | `-s 13-14` | 13, then 14A | refreshed final QC metrics and final report |
| Rebuild only the full final report/cleanup | `-s 14` | 14A | `final_results/final_result.csv` |
| Rebuild only assembly-only summary/cleanup | `--assembly-only -s 14` | 14B | `final_results/assembly_only_result.csv` |
| Add a comparison against an external assembly | add `--compare <fasta>` to any of the above | adds 14C | `final_results/compare_report/` (per-contig synteny, weak regions, alignment gaps, Circos input, optional QUAST `-R`/dnadiff/paftools/mash) |

TACO uses public step numbers **0-14**. Step 13 is final QC only. Step 14 is mode-adaptive: **14A** runs in full mode for final reporting and cleanup, while **14B** runs only when `--assembly-only` is set. **14C** runs in addition to 14A whenever `--compare` is supplied — it produces a passive contig-to-contig comparison report against `final.merged.fasta` without touching the assembly itself.

## Key Features

- Runs up to nine long-read assemblers (HiCanu, NextDenovo, Peregrine, IPA, Flye, Hifiasm, LJA, MBG, Raven) from a single command
- Supports PacBio HiFi, Oxford Nanopore, and PacBio CLR reads via `--platform`; automatically selects compatible assemblers per platform
- Input QC (Step 0): validates reads, estimates coverage, warns for low depth
- Standardizes assembly outputs for direct cross-assembler comparison
- Hybrid telomere detection with de novo k-mer discovery, built-in motif families, and per-end composite scoring
- Three-tier telomere classification: strict T2T, single-end strong, and telomere-supported
- Benchmarks all assembler outputs and the final refined assembly with BUSCO, QUAST, telomere metrics, and Merqury QV/completeness when available
- Taxon-aware scoring: BUSCO S rewarded, BUSCO D penalized, Merqury QV/completeness, telomere metrics, N50, contig count, and genome size deviation — with per-taxon weights for fungal, plant, vertebrate, insect, and other genomes
- Taxon-aware BUSCO lineage defaults: `--taxon fungal` → ascomycota, `--taxon plant` → embryophyta, etc.
- Assembly-only mode (`--assembly-only`) for publication-ready assembler benchmarking without refinement
- Conservative telomere-aware backbone refinement: D-aware duplicate filter, read-coverage diagnostic for partial overlaps, BUSCO trial validation, and "do no harm" safety comparison
- Structural quickmerge validation with parent alignment checks
- Coverage QC: sliding-window read depth analysis with GFF3 output for genome browser visualization
- Full provenance tracking: GFF3 annotation tracing every contig to its original assembler, with quickmerge region-level mapping
- Mode-aware final reporting: Step 14A creates the full final comparison report; Step 14B creates the assembly-only summary
- Optional machine-readable benchmark logs with `--benchmark`, plus decision tables and version tracking for reproducible reporting

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

**Merqury QV scoring (enabled by default):**

Merqury is automatically enabled for all platforms when `merqury.sh` and `meryl` are installed. TACO builds a reads `.meryl` k-mer database from input reads, runs Merqury on every assembler output (Step 10), and on the final refined assembly during final QC (Step 13). Per-assembly Merqury files are written under `merqury/{assembler}/` with prefixes such as `merqury/canu/canu.qv` and `merqury/canu/canu.completeness.stats`. Merqury QV and completeness are included in backbone scoring, `assemblies/assembly_info.csv`, and the final all-QC comparison table at `final_results/final_result.csv`.

By default, `--merqury-k auto` chooses an optimized k-mer size from the genome size. TACO uses Merqury's `best_k.sh` helper when it is available, otherwise it uses the same genome-size/collision-rate calculation and clamps automatic values to a practical 17-31 range for broad eukaryotic assemblies. Set `MERQURY_COLLISION_RATE` to tune the default 0.001 collision rate. Override with `--merqury-k 21` or `--merqury-k 31` when you need a fixed database for cross-run comparability.

**Accuracy note:** Merqury QV is most accurate with high-accuracy reads (PacBio HiFi or Illumina). With Nanopore or PacBio CLR reads, QV values may underestimate true assembly quality because read errors inflate k-mer error counts. Merqury completeness and relative QV ranking can still help compare assemblies generated from the same read set, but should be interpreted cautiously. TACO logs a warning when using non-HiFi reads. Disable with `--no-merqury`.

You can also provide a pre-built database:

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

## Command-Line Reference

### Parameters

| Parameter | Description |
|---|---|
| `-g`, `--genomesize` | Estimated haploid genome size (e.g., `12m`, `40m`, `2g`) |
| `-t`, `--threads` | Number of CPU threads |
| `--fastq` | Input FASTQ file (use absolute path) |
| `--taxon` | Taxonomy preset for telomere detection: `vertebrate`, `animal`, `plant`, `insect`, `fungal`, or `other` (default). Sets motif-family priors and detection behavior automatically. |
| `-m`, `--motif` | Telomere motif override (optional). Only use when the exact motif is biologically known for the species. When omitted, taxon-aware hybrid detection is used instead. |
| `--platform` | Sequencing platform: `pacbio-hifi` (default), `nanopore`, or `pacbio`. Determines compatible assemblers and the default polishing tool. |
| `-s`, `--steps` | Run selected public steps only (e.g., `1,3-5`, `13-14`). TACO uses steps 0-14; `-s 14` runs 14A unless `--assembly-only` is set. |
| `--reference`, `-ref` | Reference FASTA. Active participant: appears in QC tables, contributes its `*.telo.fasta` to the all-vs-all quickmerge candidate pool, is used as a chimera-detection alignment target during refinement, and can be selected as the backbone if it wins the scoring contest (use `--choose <assembler>` to force a specific backbone if you want to keep the reference out of selection). |
| `--compare` | Compare-only FASTA — fully passive. Goes through the same QC as the assemblers (BUSCO / Telomere / QUAST / Merqury, with a `compare` row in `final_result.csv`), and triggers **step 14C**, the contig-to-contig comparison report at `final_results/compare_report/`. The report contains: `compare_vs_final.paf` (minimap2 `-cx asm5`), `contig_to_contig.tsv` (per-compare-contig 1-to-1 mapping with identity/coverage), `unique_compare_contigs.tsv` and `unique_final_contigs.tsv` (contigs absent or < 5 % covered in the other assembly), `synteny_blocks.tsv` (1-to-1 / many-to-1 synteny summary), `weak_regions.tsv` (10 kb windows of `final.merged.fasta` < 50 % covered by compare), `compare_telomere_end_scores.tsv` and `TELOMERE_NOTES.txt` (per-contig telomere classification + threshold notes), and `circos/{karyotype.txt, links.txt, circos.conf, README.txt}` (a ready-to-render Circos plot bundle). Optional add-ons when their binaries are on PATH: `compare_quast/` (QUAST `-R compare.fa`), `compare_dnadiff/out.report` (MUMmer SNP/indel summary), `compare_paftools_variants.vcf` (`paftools.js call`), and `compare_mash_distance.tsv` (`mash dist`). Never used for backbone selection, quickmerge, telomere pool, polishing, or purge_dups. |
| `--final-fa` | Use this FASTA as the final merged assembly for steps 13/14, replacing `assemblies/final.merged.fasta`. Useful when re-running final QC on an externally produced assembly without rebuilding from step 12. |
| `--busco` | BUSCO lineage override. If omitted, TACO uses a taxon-aware default when available. |
| `--busco-download-path` | Directory where BUSCO lineage datasets are cached (passed as `--download_path` to BUSCO). Honors the `BUSCO_DOWNLOAD_PATH` env var as a fallback. |
| `--busco-offline-only` | Refuse to download BUSCO lineages over the network. If the lineage isn't already cached, BUSCO fails instead of falling back online. |
| `--choose` | Manually choose the backbone assembler |
| `--assembly-only` | Run assembler comparison only: steps 0-10, then Step 14B cleanup/reporting |
| `--auto-mode` | Backbone selection mode: `smart` (default) or `n50` |
| `--merqury` | Force-enable Merqury; if no database is provided, builds one when `meryl` is installed |
| `--merqury-db` | Enable Merqury with a specific `.meryl` database path |
| `--merqury-k` | K-mer size for an auto-built Merqury database (default: `auto`; uses Merqury `best_k.sh` when available, otherwise a genome-size/collision-rate fallback) |
| `--no-merqury` | Disable Merqury even if installed and auto-detected |
| `--benchmark` | Write optional timing/provenance files to `benchmark_logs/`; disabled by default |
| `--no-purge-dups` | Skip purge_dups after refinement |
| `--no-polish` | Skip automatic polishing after refinement |
| `--allow-t2t-replace` | Allow rescue donors to replace immutable Tier 1 (T2T) contigs. Disabled by default for safety |

### Assembly-Only Mode

Use `--assembly-only` when the goal is assembler benchmarking and comparison without refinement. TACO runs all compatible assemblers, standardizes outputs inside Step 10, runs BUSCO, telomere detection, QUAST, and optional Merqury, then writes the combined comparison table to `assemblies/assembly_info.csv` and a summary to `final_results/assembly_only_result.csv` in Step 14B.

Assembly-only mode never runs the telomere pool, refinement, or final refined-assembly QC steps. If you run `--assembly-only -s 14`, TACO runs Step 14B. If you run `-s 14` without `--assembly-only`, TACO runs Step 14A for the full final report.

Use `--benchmark` separately only when you also want machine-readable step timing and run metadata in `benchmark_logs/`.

### Benchmark Provenance Mode

`--benchmark` is for reproducibility and methods reporting, not for changing how assemblies are built. The biological benchmark is still the QC comparison itself: BUSCO, QUAST, telomere metrics, Merqury, and the final decision tables. When `--benchmark` is enabled, TACO adds a publication audit layer around that run so a future reader can reconstruct what was executed.

Use `--benchmark` when a run needs publication-ready provenance. It writes extra files in `benchmark_logs/`: exact command line, git commit/dirty state, input file size/mtime, parameters, software versions, per-step timing/status, key output file manifest, and a short methods note. FASTQ/reference SHA-256 checksums are intentionally skipped by default because large read files can be expensive to hash; set `TACO_BENCHMARK_SHA256=1` with `--benchmark` when checksums are required for an archival or paper supplement.

### Resuming And File Organization

When running selected steps with `-s`/`--steps`, TACO checks only the tools needed by those requested steps. For example, `-s 13-14` will not warn about missing assembler binaries from Steps 1-9. Before each resumed step, TACO checks for expected upstream files and prints a specific warning naming the missing file type and the earlier step that should have produced it.

Step 10 checks for raw assembler outputs from Steps 1-9 or existing normalized FASTAs. Step 12 and later check for the Step 10/11 outputs they need. Step 13 runs only final QC on the refined assembly. Step 14 does not rerun final QC; it builds the report and organizes outputs. In full mode, `-s 14` always runs 14A. In assembly-only mode, `--assembly-only -s 14` runs 14B.

For common cleanup outputs, TACO can restore active inputs from `final_results/`, `telomere_pool/`, or `temp/assemblers/` back into the working locations needed by a resumed step. As of TACO v1.3.4, restoration covers all of steps 8–14: normalized `assemblies/*.result.fasta` are restored from the corresponding `temp/assemblers/...` paths, `assembly_info.csv` and the per-merged metric CSVs are restored from `final_results/`, telomere-pool FASTAs and provenance TSVs are restored from `telomere_pool/`, and `assemblies/final.merged.fasta` is restored from `final_results/final.merged.fasta` (or from `--final-fa` when supplied — that override is authoritative). TACO v1.3.4 uses public steps 0-14; use `-s 12-14` for the full final resume path rather than older `12-17` ranges.

Cleanup keeps resumable working files in place when possible, copies stable publication-facing outputs into `final_results/`, copies telomere-pool products into `telomere_pool/`, and moves bulky transient work files into `temp/`. Final cleanup and assembly-only cleanup move raw assembler work directories into `temp/assemblers/`; normalized `assemblies/*.result.fasta` files remain the canonical comparison inputs, and Step 10 can also normalize from `temp/assemblers/` if those raw directories were already organized. If a resumed step warns that an upstream file is missing, rerun the producing step range (for example `-s 10-14`) or place the expected file back at the path shown in the warning.

## Sequencing Platform Support

TACO supports three sequencing platforms. Each assembler receives platform-appropriate flags automatically. The platform also determines the default polishing strategy: HiFi assemblies are polished with NextPolish2 by default (k-mer-based, safe for high-accuracy reads; requires `yak` for k-mer database construction), Nanopore assemblies use Medaka (neural-network polisher; falls back to Racon), and CLR assemblies use Racon.

| Platform | `--platform` | Assemblers (Steps 1-9) | Skipped |
|---|---|---|---|
| PacBio HiFi | `pacbio-hifi` (default) | canu, nextDenovo, peregrine, ipa, flye, hifiasm, lja, mbg, raven | none |
| Oxford Nanopore | `nanopore` | canu, nextDenovo, flye, raven | peregrine, ipa, hifiasm, lja, mbg |
| PacBio CLR | `pacbio` | canu, nextDenovo, peregrine, flye, raven | ipa, hifiasm, lja, mbg |

Incompatible assemblers are automatically skipped with a warning.

### Platform-Specific Assembler Compatibility

| Assembler | HiFi | Nanopore | CLR | Notes |
|---|---|---|---|---|
| Canu | `-pacbio-hifi` | `-nanopore` | `-pacbio` | All platforms supported |
| Flye | `--pacbio-hifi` | `--nano-hq` (Q20+) | `--pacbio-raw` | All platforms; set `FLYE_ONT_FLAG=--nano-raw` for older ONT |
| Hifiasm | default mode | skipped | skipped | HiFi-only as primary input |
| IPA | default mode | skipped | skipped | HiFi-only (PacBio tool) |
| Peregrine | default mode | skipped | default mode | HiFi and CLR only |
| NextDenovo | via config file | via config file | via config file | All platforms; config auto-generated |
| LJA | default mode | skipped | skipped | HiFi-only; very contiguous assemblies |
| MBG | default mode | skipped | skipped | HiFi-only; minimizer-based; included via Bioconda; uses `TACO_MBG_K` if set (default k=1501) |
| Raven | default mode | default mode | default mode | All platforms; fast OLC assembler from the `raven-assembler` package |

**Platform summary:** HiFi runs up to 9 assemblers, Nanopore runs 4 (Canu, Flye, NextDenovo, Raven), CLR runs 5 (Canu, Flye, NextDenovo, Peregrine, Raven).

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
| **Merqury** | auto (all platforms) | auto (all platforms) | auto (all platforms) | auto (all platforms) |
| **Telomere motifs** | TTAGGG + TG1-3 + Candida | TTTAGGG | TTAGGG | TTAGG (insect) / all (other) |
| **Score window** | 300 bp | 1000 bp | 1000 bp | 500 bp |
| **Backbone scoring** | S×1000 - D×600 + T2T×350 | S×1000 - D×300 + T2T×200 | S×1000 - D×500 + T2T×200 | S×1000 - D×500 + T2T×300 |
| **BUSCO trial C-drop** | 2% (strict) | 4% (relaxed) | 3% (moderate) | 2.5% (default) |
| **purge_dups mode** | two-round, upstream-default chaining | conservative single-round + polyploid warning | two-round | insect: two-round; other: conservative single-round |
| **Polishing (HiFi)** | NextPolish2 (yak k-mer based) | NextPolish2 (yak k-mer based) | NextPolish2 (yak k-mer based) | NextPolish2 (yak k-mer based) |
| **Polishing (ONT)** | Medaka → Racon | Medaka → Racon | Medaka → Racon | Medaka → Racon |
| **Polishing (CLR)** | Racon | Racon | Racon | Racon |

## Pipeline Steps

TACO exposes 15 public steps by number: **0-14**. Step 14 has two internal reporting modes, but it is still invoked as step `14` from the command line.

| Step | Description | Phase |
|---|---|---|
| 0 | Input QC and validation | Setup |
| 1-9 | Run assemblers: Canu, NextDenovo, Peregrine, IPA, Flye, Hifiasm, LJA, MBG, Raven | Assembly |
| 10 | Normalize + QC comparison: BUSCO + Telomere + QUAST + Merqury on all assemblies | QC |
| 11 | Build telomere pool (pairwise quickmerge + structural validation) | Telomere pool |
| 12 | Backbone selection and telomere-aware refinement | Refinement |
| 13 | Final QC only: BUSCO + Telomere + QUAST + Merqury on refined assembly | Final QC |
| 14 / 14A | Final comparison report + cleanup into `final_results/` (full mode) | Report |
| 14 / 14B | Assembly-only comparison summary + cleanup (only with `--assembly-only`) | Report |
| 14 / 14C | Compare-vs-final contig-to-contig report (only with `--compare`) | Report |

<img align="center" src="/docs/steps.png">

### Step 0 — Input QC

Step 0 runs automatically before assembly. It validates that the FASTQ file exists and is non-empty, parses the genome size, estimates total read bases and coverage by sampling the first 100K reads, and warns if coverage is below recommended thresholds (HiFi: 25×, ONT: 40×, CLR: 50×). It also logs which assemblers are compatible with the selected platform and confirms the BUSCO lineage setting. Step 0 runs in both full and assembly-only modes.

### Step 13 — Final QC Only

Step 13 evaluates the refined assembly produced by Step 12. It runs BUSCO, telomere detection, QUAST, and optional Merqury on `final.merged.fasta`, then writes component metric files used by Step 14A. Step 13 does not perform cleanup and does not create the full final comparison report by itself.

### Step 14 — Mode-Adaptive Reporting

Step 14 chooses its sub-mode from the run mode:

- **14A full report**: used in full refinement mode, including `-s 14`. It combines Step 10 per-assembler metrics with Step 13 final-assembly metrics, writes `final_results/final_result.csv`, and organizes final outputs.
- **14B assembly-only report**: used only when `--assembly-only` is set. It reuses Step 10 metrics, writes `final_results/assembly_only_result.csv`, and performs assembly-only cleanup.

### Full Refinement Mode (default)

Steps 0-14: runs all assemblers (1-9), normalizes and QC-compares all assemblies (10), builds the telomere pool (11), selects and refines the backbone (12), runs final QC on the refined assembly (13), and produces the final comparison report with cleanup (14A).

### Assembly-Only Mode (`--assembly-only`)

Steps 0-10, 14: runs all assemblers (1-9), normalizes and QC-compares all assemblies (10), then produces the assembly-only comparison table at `final_results/assembly_only_result.csv` (14B). Skips telomere pool (11), refinement (12), and final QC (13). This mode is designed for benchmarking and decision-making without modifying any assembly.

## Compare Mode (`--compare`)

`--compare <fasta>` adds a passive comparison against any externally produced assembly (e.g. an NCBI GenBank genome of the same or a closely related strain). The compare FASTA is **never** used for backbone selection, telomere-pool quickmerge, polishing, or `purge_dups`. It is treated as a read-only reference whose only effect on TACO's run is to land in the QC tables and trigger an extra reporting sub-step.

`--compare` works in two invocation modes and produces identical outputs in both:

- **Full pipeline:** `taco ... --compare X.fa` runs steps 0–14 normally; the comparison report is generated as part of step 14.
- **Resume on an existing run:** `taco ... --compare X.fa -s 10,13,14` after a prior full run reuses the existing assemblies and merged FASTA and only adds the comparison artefacts.

### What runs where

| Stage | Effect of `--compare` |
|---|---|
| Step 10 (Normalize + QC) | Compare FASTA is normalized to `assemblies/compare.result.fasta` and gets BUSCO, telomere, QUAST, and Merqury rows in `assembly_info.csv` and (later) `final_result.csv`. |
| Steps 11 / 12 | **Skipped for compare** — it is excluded from the telomere pool, the all-vs-all quickmerge candidate set, polish, and purge_dups. |
| Step 13 (Final QC) | Unchanged — runs on `final.merged.fasta` only. |
| **Step 14C** (new) | Builds the contig-to-contig report at `final_results/compare_report/` against `final.merged.fasta`, then proceeds to standard 14A cleanup. |

### Outputs in `final_results/compare_report/`

Step 14C aligns the compare FASTA to `final.merged.fasta` with `minimap2 -cx asm5 --secondary=no` and writes the following layered diagnostics:

| File | Contents |
|---|---|
| `compare_vs_final.paf` | Raw minimap2 alignment (one row per alignment block). All downstream files are derived from this. |
| `contig_to_contig.tsv` | One row per compare contig: `best_target_contig`, `best_target_aligned_bases`, `best_target_coverage_pct`, total `aligned_bases`, average `identity_pct`, `n_significant_targets`, `relationship` (`1-to-1` / `1-to-N (split)`), `all_targets_breakdown` (every significant target with its bp and coverage), and per-end telomere flags for both the compare contig and its best target. |
| `contig_to_contig_pairs.tsv` | One row per (compare_contig, final_contig) pair above the significance threshold (max of 20 % of the compare contig OR 100 kb absolute). Per-pair compare coverage, final coverage, identity, dominant strand, four boundary-touch flags, and four pair-localized telomere flags. This is what makes 1-to-N splits readable. |
| `synteny_blocks.tsv` | At-a-glance synteny summary derived from the per-pair table — labels each block as `1-to-1`, `1-to-many (compare splits across N final contigs)`, `many-to-1 (M compare contigs → 1 final)`, or `many-to-many`. |
| `split_mappings.tsv` | Only compare contigs that map to >1 significant final contig. Includes a `chimera_evidence` text column that classifies the split using the telomere pattern at outer + inner boundaries: full-coverage split, strong chimera signal, chimeric extension into a fragment, or split-likely-real. |
| `unique_compare_contigs.tsv` | Compare contigs with < 5 % coverage in `final.merged.fasta` (or no alignment). |
| `unique_final_contigs.tsv` | Mirror file: `final.merged.fasta` contigs with < 5 % coverage from `--compare`. |
| `weak_regions_final.tsv` / `.gff3` | 10 kb windows of `final.merged.fasta` < 50 % covered by compare. (`weak_regions.tsv` is kept as a stable alias for back-compat.) |
| `weak_regions_compare.tsv` / `.gff3` | 10 kb windows of `--compare` < 50 % covered by `final.merged.fasta`. Empty for high-coverage runs — see the next file for fine-grained gap detection. |
| `final_alignment_gaps.tsv` / `.gff3` | Every uncovered interval ≥ 100 bp on `final.merged.fasta`, with `kind` ∈ {`leading`, `trailing`, `internal`, `internal_junction_candidate`, `*_tip`, `full_contig_unaligned`}. Surfaces telomere-padding tips, chimera-broken tails, and TACO-unique sequence at sub-window resolution. |
| `compare_alignment_gaps.tsv` / `.gff3` | Mirror file on the compare side. The chimera-junction bridge inside a 1-to-many compare contig appears here as a small `internal_junction_candidate` (typically 100–500 bp). |
| `compare_telomere_end_scores.tsv` | Copy of `assemblies/compare.telomere_end_scores.tsv` from step 9, surfaced inside the report for convenience. |
| `TELOMERE_NOTES.txt` | Plain-text summary of compare telomere classifications + thresholds + notes on widening `--telo-score-window` if telomere repeats are detected just inside the contig ends. |
| `circos/` | Self-contained Circos input bundle: `karyotype.txt` (compare = `C_<name>`, final = `F_<name>`, distinct palettes), `links.txt` (one ribbon per minimap2 block ≥ 5 kb), `circos.conf`, and `README.txt` with run instructions. Render with `cd circos && circos -conf circos.conf`. |
| `compare_quast/` | *Optional* — QUAST run with `-R compare.fa` against `final.merged.fasta` for genome fraction, NA50, and reference-based misassembly counts. Skipped silently when QUAST is not on PATH. |
| `compare_dnadiff/` | *Optional* — MUMmer `dnadiff` SNP/indel summary. Skipped silently when `dnadiff` is not on PATH (install MUMmer4 to enable). |
| `compare_paftools_variants.vcf` | *Optional* — `paftools.js call -L 10000` over the sorted PAF for assembly-vs-assembly SNVs and SVs. Ships with `minimap2`; skipped silently when `paftools.js` / `k8` aren't on PATH. |
| `compare_mash_distance.tsv` | *Optional* — `mash dist compare.fa final.fa` scalar distance. Skipped silently when `mash` is not on PATH. |

### How to read the report

The report is layered so that each file answers a more specific question than the one before it:

1. Start with **`split_mappings.tsv`**. If it has any rows, those are the chimera or split candidates worth investigating. The `chimera_evidence` column tells you which kind.
2. For each split, open **`contig_to_contig_pairs.tsv`** and filter on the compare contig name. The four boundary-touch flags + four telomere-at-pair flags reveal whether each pair sits on compare's outer end or near the join, and whether the corresponding final-contig boundary carries a telomere.
3. **`final_alignment_gaps.tsv`** explains every region of TACO's assembly that the compare doesn't cover — telomere arrays the compare lost, chimera-broken chromosome ends, and minor alignment-edge artifacts. **`compare_alignment_gaps.tsv`** does the mirror, and crucially surfaces the small (100–500 bp) chimera-junction bridges that the window-based `weak_regions_compare.tsv` misses.
4. For visual inspection, render `circos/` with Circos, or load `weak_regions_*.gff3` and `*_alignment_gaps.gff3` as IGV / JBrowse tracks alongside the FASTA. A chimera looks like two ribbons from the same `C_<name>` ideogram converging on different `F_<name>` ideograms with opposing strands.

### `--final-fa` companion flag

If you want to run step 14C against an externally-produced final assembly without rebuilding from step 12, pass `--final-fa <path>`. TACO copies that FASTA into `assemblies/final.merged.fasta` (overriding any cached copy) and runs steps 13/14 against it. Combined with `--compare`, this lets you compare any two assemblies head-to-head without re-running the full pipeline:

```bash
taco -g 40m -t 30 --fastq reads.fastq --platform nanopore --taxon fungal \
  --final-fa /path/to/your_final.fasta \
  --compare  /path/to/external.fasta \
  -s 13,14
```

### Recommended invocations

```bash
# Full pipeline + comparison against an NCBI reference
taco -g 40m -t 30 --fastq reads.fastq --platform nanopore --taxon fungal \
  --busco-download-path /shared/busco_downloads --busco-offline-only \
  --compare /path/to/GCA_xxxxx.fna

# Add a comparison after a previous full run (no reassembly)
cd <previous_run_dir>
taco -g 40m -t 30 --fastq reads.fastq --platform nanopore --taxon fungal \
  --compare /path/to/GCA_xxxxx.fna \
  -s 10,13,14

# Compare two existing assemblies without rerunning
taco -g 40m -t 30 --fastq reads.fastq --platform nanopore --taxon fungal \
  --final-fa /path/to/your.fasta \
  --compare  /path/to/other.fasta \
  -s 13,14
```

## Telomere Detection

TACO v1.3.4 uses a taxon-aware hybrid telomere detection system that combines built-in motif families with de novo k-mer discovery.

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

TACO ranks assemblies using a composite score. Merqury QV and completeness are included by default when `merqury.sh` + `meryl` are installed (otherwise these terms contribute 0):

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

A D-rise (duplicated BUSCO increase) check catches cases where a rescue introduces redundant copies of single-copy orthologs — a sign of retained haplotigs or mis-joined contigs. All thresholds can be overridden via `STEP12_MAX_BUSCO_C_DROP`, `STEP12_MAX_BUSCO_M_RISE`, and `STEP12_MAX_BUSCO_D_RISE` environment variables. Legacy `STEP13_*` names are still accepted for compatibility with older run scripts.

Additional environment variables for fine-tuning Step 12: `STEP12_MAX_ACCEPTED`, `STEP12_MIN_BP_RATIO`, `STEP12_BUSCO_TRIAL_TIMEOUT` (seconds per BUSCO trial attempt, default 43200; set 0 to disable), `STEP12_BUSCO_ALLOW_DOWNLOAD` (default 0; set 1 to permit online BUSCO lineage fallback), `STEP12_SKIP_BUSCO_TRIAL` (set 1 to use structural rescue checks only), `CHIMERA_MIN_CROSS_COV` (minimum cross-assembly coverage for chimera mapping check, default 0.60), `SELFDEDUP_ENABLE` plus `SELFDEDUP_COV` / `SELFDEDUP_ID` (optional self-dedup), `RESCUE_MIN_IDENT`, `RESCUE_MIN_ALN_BP`, `RESCUE_MIN_COV_BB`, `RESCUE_MIN_COV_DONOR`, `RESCUE_MIN_EXT`, `NOVEL_DUP_COV`, `NOVEL_DUP_ID`, `NOVEL_UPGRADE_TCOV`, `NOVEL_MAX_D_RISE`, and `PARTIAL_T2T_MIN_TCOV` / `PARTIAL_T2T_MIN_BP` / `PARTIAL_T2T_MIN_QCOV` / `PARTIAL_T2T_MIN_IDENT` for coverage-guided partial T2T replacement.

### Taxon-Specific Rescue Limits (Step 12F)

The maximum number of accepted rescue candidates per run is taxon-aware: fungi allow up to 20 rescues (many small chromosomes), vertebrates 10, plants 8 (conservative due to polyploidy risk), and other taxa 15. This prevents runaway replacement in complex genomes.

### D-Aware Duplicate Filter (Step 12F2)

Before adding "novel" pool T2T contigs, TACO checks each candidate against the current backbone with minimap2 and applies a five-tier decision:

1. **Overlaps Tier 1 (T2T) backbone** → reject (pure duplicate, backbone already has T2T).
2. **Overlaps Tier 2 (non-T2T) backbone at ≥80% target coverage** → upgrade (replace backbone with T2T — better telomere evidence, comparable size).
3. **Overlaps Tier 2 at 50–80% target coverage** → read-coverage diagnostic. TACO maps the input reads to the backbone contig with a platform-specific minimap2 preset and compares median coverage in the T2T-covered region vs the uncovered region. If the uncovered region has very low coverage (< 30% of covered), the backbone is chimeric and the T2T is the real chromosome — TACO replaces the backbone. If coverage is normal, the backbone is real — TACO rejects the novel contig (adding it would increase BUSCO D). Configurable via `CHIMERIC_COV_RATIO` (default 0.30).
4. **Overlaps Tier 2 at < 50%** → reject (insufficient overlap for any useful decision).
5. **No significant overlap** → add as genuinely novel chromosomal region (with optional BUSCO D check: `NOVEL_MAX_D_RISE`).

Thresholds: `NOVEL_DUP_COV` (default 0.80), `NOVEL_DUP_ID` (default 0.90), `NOVEL_UPGRADE_TCOV` (default 0.80).

### Taxon-Specific purge_dups Behaviour (Step 12H)

purge_dups strategy is taxon-aware. Fungi use two-round purging (`-2`) with upstream-default chaining lengths, which is safer than short-match fungal tuning when the selected backbone is already close to the expected genome size. Vertebrate, animal, and insect genomes use two-round purging for high-heterozygosity haplotig/overlap cleanup. Plant genomes use conservative single-round purging with stricter sequence-level evidence to avoid collapsing homeologous sequences in polyploid species — use `--no-purge-dups` or `PURGE_DUPS_MODE=skip` if this is still too aggressive. TACO runs `get_seqs -e` by default, so only end duplications are removed, and rejects purged output if the taxon-specific size drop or expected genome-size floor suggests over-purging. If purging makes an overlarge assembly substantially closer to the expected genome size without falling below the floor, TACO accepts that drop instead of preserving duplicated sequence. Coverage cutoffs from `calcuts` are logged for debugging; override with `PURGE_DUPS_CALCUTS`. Advanced overrides: `PURGE_DUPS_MODE` (`auto`, `single`, `two-round`, `skip`), `PURGE_DUPS_EXTRA_OPTS`, `PURGE_DUPS_GET_SEQS_OPTS`, `PURGE_DUPS_MAX_BP_DROP`, and `PURGE_DUPS_MIN_EXPECTED_RATIO`.

### N50-only Mode

`--auto-mode n50` selects the assembly with the highest N50. This reproduces legacy behavior but may favor contiguous assemblies that lack completeness.

## Step 12 — Backbone-First Telomere-Aware Refinement

Step 12 (backbone refinement) adopts a backbone-first assembly philosophy with a **two-tier confidence model**:

- **Tier 1 (Immutable)**: T2T contigs — contigs with verified telomere signal at both ends. These are treated as protected chromosomal anchors and are never replaced during rescue, unless `--allow-t2t-replace` is explicitly set. This protects the highest-confidence contigs from accidental degradation.
- **Tier 2 (Editable)**: Backbone contigs — gap-fill contigs that cover chromosomal regions not represented by T2T contigs. These may be replaced by telomere-bearing rescue donors if the replacement passes BUSCO trial validation.

Backbone contigs are preserved by default to maintain BUSCO completeness — purge_dups at Step 12H handles haplotig removal using read-coverage evidence. Novel T2T additions undergo a D-aware duplicate filter that either upgrades Tier 2 backbone contigs (replacing non-T2T with T2T) or rejects pure duplicates of Tier 1 contigs. Rescue donors must carry verified telomere signal.

### Step 12 Sub-step Flow

1. **12A** — refresh Merqury QV scoring if needed.
2. **12B** — auto-select backbone assembler from the Step 10 comparison table.
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
11. **12K** — final assembly coverage QC: maps the input reads to the final assembly with a platform-specific minimap2 preset, computes sliding-window coverage (default 5 kb), and flags zero-coverage gaps, very-low-coverage regions, and sudden coverage drops. Outputs `coverage_summary.tsv`, `weak_regions.tsv`, and `weak_regions.gff3` (loadable in IGV).
12. **12L** — "do no harm" safety comparison: compares final assembly vs original backbone for size, telomere count, and genome size deviation. If refinement degraded quality, both assemblies are preserved with a `refinement_warning.txt` explaining the issues.

### BUSCO Trial Validation

TACO validates each telomere rescue candidate by building a trial assembly where one backbone contig is replaced by one donor contig, then running BUSCO with the same lineage selected by the user. Rejection is triggered by three independent BUSCO metrics: C% drop (completeness loss), M% rise (missing gene increase), and D% rise (duplicated BUSCO increase, catching retained haplotigs). Rejection thresholds are taxon-aware: fungi use strict thresholds (2% C-drop, 2% D-rise), plants use relaxed thresholds (4% C-drop, 6% D-rise, accounting for polyploidy), and vertebrates use moderate thresholds (3% C-drop, 4% D-rise). An additional safety check rejects `replace_single_with_better` candidates if telomere evidence weakens at either end after replacement (suspicious size drops >30% also trigger rejection). This greedy, sequential approach ensures that each accepted rescue improves or maintains assembly quality. The trial summary TSV includes `replacement_class` and `D` (duplicated %) columns for full traceability.

Trial BUSCO stdout/stderr logs are written to `assemblies/rescue_trials/*.busco.*.log`; final BUSCO logs are written to `busco/*.busco.*.log`. Each BUSCO attempt has a timeout controlled by `STEP12_BUSCO_TRIAL_TIMEOUT` (default 43200 seconds; set to 0 to disable). TACO tries BUSCO validation in offline mode by default and does not fall back to online lineage download/checks unless `STEP12_BUSCO_ALLOW_DOWNLOAD=1`, which avoids cluster runs hanging in BUSCO startup. If BUSCO is unavailable or wedged, set `STEP12_SKIP_BUSCO_TRIAL=1` to resume Step 12 with structural rescue checks only.

### Post-Refinement Stack

After the rescued/combined assembly is produced, TACO runs purge_dups by default to remove leftover haplotigs, overlapping fragments, and residual duplicates. purge_dups strategy is taxon-aware and guarded against over-purging: fungi use two-round purging with upstream-default chaining lengths; vertebrate, animal, and insect genomes use two-round purging; plants use conservative single-round purging to preserve homeologs. TACO keeps the unpurged assembly if `dups.bed` is empty or if the purged output fails the taxon-specific size safety check. Coverage cutoffs from `calcuts` are logged; override with `PURGE_DUPS_CALCUTS` if automatic thresholds are incorrect. This is followed by automatic polishing selected from `--platform`: HiFi assemblies use NextPolish2 (maps HiFi reads to assembly with minimap2, sorts BAM with samtools, builds yak k-mer databases, then runs k-mer-based correction), Nanopore uses Medaka (falls back to Racon), and CLR uses Racon. Both steps can be skipped with `--no-purge-dups` and `--no-polish`.

Strict T2T pool contigs only replace backbone contigs directly when they cover most of the target backbone contig. For smaller but substantial partial T2T hits, TACO runs a read-coverage diagnostic on the unmatched backbone sequence. If the unmatched region has low coverage, the backbone is treated as duplicated/chimeric and replaced by the T2T contig; if coverage is normal or inconclusive, TACO keeps the backbone to preserve BUSCO gene content. Tiny telomere-repeat fragments are rejected as partial replacement donors.

### Diploid and Polyploid Note

TACO is designed for producing a best primary-style chromosome-level assembly, not a fully phased diploid or polyploid reconstruction. For strongly diploid or polyploid genomes, telomere-bearing contigs from different haplotypes may appear as rescue donors, and purge_dups may collapse alternative haplotigs. This is acceptable when the goal is a cleaned primary reference assembly.

### Provenance GFF3

TACO writes a GFF3 annotation file (`final.merged.provenance.gff3`) alongside the final assembly. Each contig gets one GFF3 record (type=contig) spanning its full length, with attributes documenting its full provenance chain: `source_assembler` (which assembler produced it), `assembler_contig` (the original contig name from that assembler before Step 11 pool renaming), `source_type` (assembler or quickmerge), `role` (backbone, upgrade_donor, or novel_t2t), `replacement_class` (for upgrade donors), `replaced_contig` (which backbone contig was replaced), and `description` (a human-readable summary like "Entire replacement: peregrine contig 'contig_5' replaced by canu contig 'tig00000015' (class: upgrade_tier2_to_t2t)").

For quickmerge-derived contigs, the GFF3 includes additional contig-level attributes (`qm_assembler1`, `qm_assembler2`) identifying the two source assemblers, plus child records (type=region) with `Parent` linking to the contig. Each region record spans the genomic coordinates contributed by a specific assembler, with `source_assembler` and `assembler_contig` showing the original source. For example, a quickmerge contig produced from canu × flye will have region records like "Region 1-500000 from canu contig 'tig00001'" and "Region 400000-900000 from flye contig 'contig_3'", enabling users to trace every base pair back to its assembler of origin.

A companion file `pool_contig_provenance.tsv` maps every pool contig back to its source assembler and original contig name, with extended columns for quickmerge contigs: `qm_assembler1`, `qm_assembler2`, and `qm_regions` (semicolon-delimited `start-end:assembler:contig` entries).

### Coverage QC GFF3 and Reports

TACO maps the input reads back to the final assembly (Step 12K) and scans for assembly errors using a platform-specific minimap2 preset and a sliding-window coverage analysis (default 5 kb window). Three output files are produced:

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

## Outputs And Reports

### Key Output Files

| File | Produced by | Purpose |
|---|---|---|
| `assemblies/assembly_info.csv` | Step 10 | Unified assembler comparison table used for benchmarking and backbone selection |
| `final_results/assembly_only_result.csv` | Step 14B | Publication-facing summary for `--assembly-only` runs |
| `final_results/final_assembly.fasta` | Step 14A | Final refined candidate assembly for downstream analysis |
| `final_results/final_result.csv` | Step 14A | Final report combining assembler metrics with the final refined assembly metrics |
| `final_results/final.merged.provenance.gff3` | Step 12/14A | Contig-level and region-level provenance for the refined assembly |
| `final_results/coverage_summary.tsv` | Step 12K/14A | Per-contig read-depth QC summary for the final assembly |
| `final_results/weak_regions.tsv` and `.gff3` | Step 12K/14A | Coverage-warning regions for manual inspection in IGV/JBrowse |
| `benchmark_logs/` | `--benchmark` | Optional reproducibility metadata, software versions, timing, and manifest files |

```
project_directory/
├── assemblies/
│   ├── assembly_info.csv                # Unified comparison table
│   ├── canu.result.fasta                # Normalized assembly outputs
│   ├── nextDenovo.result.fasta
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
│   ├── final_result.csv                 # Full-mode final comparison report (Step 14A)
│   ├── final.merged.fasta               # Refined assembly used internally for final QC
│   ├── final_assembly.fasta             # Publication-facing refined assembly copy
│   ├── final.merged.provenance.gff3     # GFF3 provenance: full assembler tracing per contig
│   ├── pool_contig_provenance.tsv       # Pool contig → assembler + original name mapping
│   ├── quickmerge_validation.tsv        # Quickmerge structural validation decisions
│   ├── telomere_pool_decisions.tsv      # Per-contig pool classification decisions
│   ├── rescue_trial_summary.tsv         # BUSCO rescue trial results
│   ├── rescue_rejection_summary.txt     # Rescue rejection summary
│   ├── single_tel.candidates.tsv        # Candidate single-telomere rescue contigs
│   ├── selection_decision.txt           # Backbone-selection formula and selected assembler
│   ├── coverage_summary.tsv             # Per-contig coverage stats (median, mean, zero/low bp)
│   ├── weak_regions.tsv                 # Flagged weak windows with coords, source assembler, flag
│   ├── weak_regions.gff3                # GFF3 coverage warnings: load in IGV to see weak spots
│   ├── {backbone}.backbone.original.fasta  # Original backbone (for do-no-harm comparison)
│   ├── refinement_warning.txt           # Quality warnings if refinement degraded backbone
│   └── assembly_only_result.csv         # Assembly-only comparison summary (Step 14B)
├── telomere_pool/                       # Structured copy of telomere-pool intermediates
│   ├── protected_telomere_contigs.fasta
│   ├── t2t_clean.fasta
│   ├── single_tel_best_clean.fasta
│   ├── telomere_supported_best_clean.fasta
│   ├── pool_contig_provenance.tsv
│   └── telomere_pool_decisions.tsv
├── quast_out/                           # Pre-refinement QUAST output
├── quast_final/                         # Final-assembly QUAST output
├── logs/                                # Per-step log files
├── temp/                                # Organized temporary alignment/work files after cleanup
│   ├── merge/
│   ├── assemblers/                      # Raw assembler work directories after cleanup
│   ├── telomere/
│   ├── polish/
│   ├── purge_dups/
│   ├── qc/
│   ├── busco/
│   └── log/
├── benchmark_logs/                      # Optional timing/provenance data, only with --benchmark
│   ├── run_metadata.tsv
│   ├── run_manifest.json
│   ├── software_versions.tsv
│   ├── step_benchmark.tsv
│   ├── output_manifest.tsv
│   ├── methods_note.txt
│   └── run_summary.txt
└── version.txt                          # Software versions
```

## Project Structure

```
TACO/
├── setup.py                # pip install entry point
├── run_taco                # Shell wrapper (no install needed)
├── taco/                   # Python package
│   ├── __init__.py         # Package metadata (v1.3.4)
│   ├── __main__.py         # CLI entry point: taco [options]
│   ├── cli.py              # Argument parsing
│   ├── pipeline.py         # Pipeline runner, logging, benchmarking
│   ├── steps.py            # Step implementations (0-14, including 14A/14B)
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

**Raven skipped even after installing `raven-assembler`:** Confirm that the environment where TACO is running is the same environment where Raven was installed. Activate the environment and check `which raven`. With micromamba, use `micromamba activate taco` before running TACO, or make sure `/path/to/env/bin` is on `PATH`. The Bioconda package is `raven-assembler`, but the executable TACO looks for is `raven`.

**Telomere motif appears incorrect:** Do not force `--motif` unless the telomere repeat is biologically known for your species. Use `--taxon` to select the appropriate preset instead. TACO's built-in motif families cover canonical TTAGGG, budding yeast TG1-3, Candida, plant TTTAGGG, and insect TTAGG repeats.

**purge_dups or polishing not running:** These tools must be installed in the conda environment. Use `conda install -c bioconda purge_dups nextpolish2 yak racon medaka` or skip with `--no-purge-dups` / `--no-polish`. For HiFi polishing, NextPolish2 and yak are both required. For Nanopore polishing, Medaka is preferred; if unavailable, Racon is used as fallback.

**`taco: command not found`:** Activate the TACO conda/micromamba environment and reinstall the package with `pip install -e .` from the repository root. If you prefer the wrapper, run `./run_taco` from the repository directory.

**Missing Python modules:** TACO uses only the Python standard library. If you see import errors, ensure Python >= 3.8 is installed and that you are running the installed `taco` command or the repository-local `./run_taco` wrapper.

**Merqury not working:** Merqury is enabled by default when `merqury.sh` + `meryl` are installed. The reads `.meryl` database is built automatically from input reads. If you have a pre-built database, use `--merqury-db path/to/reads.meryl`. Install with `conda install -c bioconda merqury meryl`. TACO writes organized outputs such as `merqury/canu/canu.qv` and also searches legacy flat prefixes such as `merqury/canu.qv`; if `assembly.merqury.csv` is empty, check the step log for the warning that lists which Merqury files were found. Nanopore and PacBio CLR runs log a warning because QV from non-high-accuracy reads can be underestimated; completeness and relative assembler ranking are still reported but should be interpreted cautiously. Disable with `--no-merqury`.

**Merqury completeness exists but QV is `NA`:** This means Merqury produced a completeness file but no parseable `.qv` file for that assembly. TACO reports `NA` rather than leaving the final report blank. Check the corresponding `merqury/{label}/` directory and step log to confirm whether Merqury wrote a QV table under a non-standard name.

**`-s 14` produced the full report instead of assembly-only output:** This is expected. Step 14A is the default full-mode report. Step 14B runs only with `--assembly-only`, for example `taco ... --assembly-only -s 14`.

## Citation

If you use TACO in a publication, please cite the software and archive the exact release used for reproducibility (e.g., via Zenodo).

TACO was developed at the Grainger Bioinformatics Center, Field Museum of Natural History, Chicago, Illinois, USA.

## Changelog

See [docs/CHANGELOG.md](docs/CHANGELOG.md) for the full version history, including detailed notes on the v0.5.6 → v1.0.0 Bash-to-Python conversion.

## License

TACO is released under the [MIT License](LICENSE).
