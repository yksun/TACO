<img align="center" src="/docs/taco-icon.png">

# TACO

**Telomere-Aware Contig Optimization**

TACO is a telomere-aware all-in-one multi-assembler comparison and refinement pipeline for genome assembly benchmarking, decision-making, and chromosome-end improvement. Developed for small eukaryotic genomes with a focus on fungal genomes, TACO runs multiple assemblers, standardizes their outputs, evaluates assembly quality, detects telomere-supported contigs, and can either **(1)** stop after generating a unified assembler comparison table for benchmarking, or **(2)** continue into telomere-aware backbone refinement for an improved chromosome-scale candidate assembly.

TACO was developed at the **Grainger Bioinformatics Center, Field Museum of Natural History**.

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
  -m TTAGGG
```

**Assembly-only comparison** (benchmarking only):

```bash
taco -g 12m -t 16 \
  --fastq /path/to/reads.fastq \
  -m TTAGGG \
  --assembly-only
```

**With Merqury QV scoring:**

```bash
taco -g 12m -t 16 \
  --fastq /path/to/reads.fastq \
  -m TTAGGG \
  --merqury-db reads.meryl
```

## Usage

### Parameters

| Parameter | Description |
|---|---|
| `-g`, `--genomesize` | Estimated haploid genome size (e.g., `12m`, `40m`, `2g`) |
| `-t`, `--threads` | Number of CPU threads |
| `--fastq` | Input FASTQ file (use absolute path) |
| `-m`, `--motif` | Telomere motif or seed motif |
| `--platform` | Sequencing platform: `pacbio-hifi` (default), `nanopore`, or `pacbio` |
| `-s`, `--steps` | Run selected steps only (e.g., `1,3-5`) |
| `--fasta` | External pre-assembled FASTA to include in comparison. Also used as the reference genome for Redundans reference-guided scaffolding in Step 12. When omitted, Step 12 runs in pure de novo mode (long-read-only scaffolding). |
| `--busco` | Run BUSCO (optionally specify lineage dataset) |
| `--choose` | Manually choose the backbone assembler |
| `--assembly-only` | Stop after assembler comparison |
| `--auto-mode` | Backbone selection mode: `smart` (default) or `n50` |
| `--merqury` | Enable Merqury with auto-detected `.meryl` database |
| `--merqury-db` | Enable Merqury with a specific `.meryl` database path |
| `--no-merqury` | Disable Merqury even if available |

### Assembly-Only Mode

Use `--assembly-only` when the goal is assembler benchmarking and comparison without refinement. TACO runs all assemblers, standardizes outputs, runs BUSCO, telomere detection, QUAST, and optional Merqury, then writes the combined comparison table to `assemblies/assembly_info.csv` and a summary to `final_results/assembly_only_result.csv`.

## Sequencing Platform Support

TACO v1.0.0 supports three sequencing platforms. Each assembler receives platform-appropriate flags automatically.

| Platform | `--platform` | Assemblers Used | Notes |
|---|---|---|---|
| PacBio HiFi | `pacbio-hifi` (default) | canu, nextDenovo, peregrine, IPA, flye, hifiasm | All 6 assemblers |
| Oxford Nanopore | `nanopore` | canu, nextDenovo, flye, hifiasm | Peregrine and IPA skipped |
| PacBio CLR | `pacbio` | canu, nextDenovo, peregrine, flye, hifiasm | IPA skipped |

Incompatible assemblers are automatically skipped with a warning. For example, IPA only supports PacBio HiFi reads, and Peregrine does not support Nanopore reads.

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

With `--assembly-only`, TACO follows the comparison path (Steps 1-11, 18) and stops before backbone refinement.

## Telomere Detection

TACO v1.0.0 uses a hybrid telomere detection system that combines built-in motif families with optional de novo k-mer discovery.

### Built-in Motif Families

TACO ships with three motif families covering the most common telomere repeats in eukaryotes: the canonical vertebrate/filamentous fungal TTAGGG repeat, the budding yeast TG1-3/C1-3A degenerate repeat, and the Candida 23-bp repeat (ACGGATGTCTAACTTCTTGGTGT). Users can specify a motif with `-m` or allow automatic detection.

### Scoring System

Each contig end is scored using a composite of four metrics: telomere density (weight 0.40), longest consecutive run (weight 0.30), distance of repeats from the contig terminus (weight 0.20), and covered base pairs (weight 0.10). The composite score ranges from 0 to 1.

### Classification Tiers

Contigs are classified into three tiers based on their end scores: **strict T2T** contigs have strong telomere signal at both ends (score >= 0.25 at each end), **single-end strong** contigs have strong signal at one end only, and **telomere-supported** contigs have at least weak signal (score >= 0.08) at one end.

## Assembly Selection Strategy

When `--choose` is not provided, TACO automatically selects the backbone assembly for refinement.

### Smart Scoring (default)

TACO ranks assemblies using a composite score that prioritizes biological completeness and chromosome-end support:

```
score = BUSCO_S * 1000 + T2T * 300 + single_tel * 150
      + MerquryComp * 200 + MerquryQV * 20
      - contigs * 30 + log10(N50) * 150
```

BUSCO single-copy completeness (S%) is used instead of total completeness (C%) to avoid rewarding highly duplicated assemblies. Telomere metrics contribute meaningfully but do not dominate, since the telomere pool can rescue missing ends during refinement. Merqury metrics are included when available.

### N50-only Mode

`--auto-mode n50` selects the assembly with the highest N50. This reproduces legacy behavior but may favor contiguous assemblies that lack completeness.

## Step 12 — Backbone Refinement with Redundans

Step 12 refines the backbone assembly in several sub-steps: redundancy reduction against the protected T2T pool, Redundans-based reduction of remaining heterozygous duplicates, telomere rescue, and Redundans scaffolding + gap closing with the original long reads.

### Redundans Integration

TACO uses [Redundans](https://github.com/Gabaldonlab/redundans) (Pryszcz & Gabaldón 2016) for three purposes in Step 12:

1. **Redundancy reduction** (Step 12D Pass 2) — After removing backbone contigs near-identical to T2T contigs (95%/95%), Redundans detects and removes heterozygous/duplicate contigs among the surviving backbone. Thresholds: identity >= 0.51, overlap >= 0.80 (configurable via `RED_IDENTITY` / `RED_OVERLAP` environment variables).

2. **Long-read scaffolding** (Step 12G2) — The combined assembly (protected T2T + backbone) is scaffolded using the original sequencing reads (HiFi, ONT, or CLR). This can join fragmented backbone contigs where read evidence supports a join.

3. **Gap closing** (Step 12G2) — Gaps introduced during scaffolding are filled using the same long reads.

### Reference-Guided vs De Novo Mode

The `--fasta` parameter controls whether Redundans operates in reference-guided or de novo mode:

- **With `--fasta`**: The external FASTA is treated as a reference genome and passed to Redundans via `-r`. This enables reference-guided scaffolding where Redundans uses the reference chromosome structure to order and orient contigs, in addition to the long-read evidence.

- **Without `--fasta`** (pure de novo): Scaffolding relies solely on long-read evidence to join contigs. No external reference is used. This is the standard mode for de novo genome assembly projects.

If `redundans.py` is not installed, TACO falls back to a minimap2-based fragment removal for reduction and skips scaffolding/gap closing.

## Output Structure

```
project_directory/
├── assemblies/
│   ├── assembly_info.csv          # Unified comparison table
│   ├── canu.fasta                 # Normalized assembly outputs
│   ├── nextdenovo.fasta
│   ├── ...
│   └── *.busco/                   # BUSCO results per assembly
├── final_results/
│   ├── final_result.csv           # Final comparison report
│   ├── final_assembly.fasta       # Refined assembly (full mode)
│   └── assembly_only_result.csv   # Comparison summary (assembly-only)
├── telomere_pool/                 # Telomere pool intermediates
├── quast_results/                 # QUAST output
├── logs/                          # Per-step log files
├── benchmark_logs/                # Machine-readable benchmark data
│   ├── run_metadata.tsv
│   ├── step_benchmark.tsv
│   └── run_summary.txt
└── version.txt                    # Software versions
```

## Project Structure

```
TACO/
├── setup.py                # pip install entry point
├── run_taco                # Shell wrapper (no install needed)
├── taco/                   # Python package
│   ├── __init__.py         # Package metadata (v1.0.0)
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

**Telomere motif appears incorrect:** Do not assume the same motif for all fungi. Choose a motif appropriate for the target organism. TACO's built-in motif families cover canonical TTAGGG, budding yeast TG1-3, and Candida repeats.

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
