<img align="center" src="/docs/taco-icon.png">

# TACO

**Telomere-Aware Contig Optimization**

TACO is a telomere-aware all-in-one multi-assembler comparison and refinement pipeline for genome assembly benchmarking, decision-making, and chromosome-end improvement. Developed for small eukaryotic genomes with a focus on fungal genomes, TACO runs multiple assemblers, standardizes their outputs, evaluates assembly quality, detects telomere-supported contigs, and can either **(1)** stop after generating a unified assembler comparison table for benchmarking, or **(2)** continue into telomere-aware backbone refinement for an improved chromosome-scale candidate assembly.

TACO was developed at the **Grainger Bioinformatics Center, Field Museum of Natural History**.

![Latest Version](https://img.shields.io/github/v/tag/ysun-fieldmuseum/TACO?label=Latest%20Version)
![Last Commit](https://img.shields.io/github/last-commit/ysun-fieldmuseum/TACO)
![Issues](https://img.shields.io/github/issues/ysun-fieldmuseum/TACO)
![License](https://img.shields.io/github/license/ysun-fieldmuseum/TACO)

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

# Run TACO
python3 -m taco --help
```

### Requirements

TACO requires a Unix-like system (Linux or macOS) with Python >= 3.8 and Conda. All pipeline dependencies (assemblers, analysis tools) are specified in `taco-env.yml`. TACO's Python modules use only the standard library with no additional pip packages.

## Quick Start

**Full refinement run** (comparison + telomere-aware refinement):

```bash
mkdir -p my_project && cd my_project

python3 -m taco -g 12m -t 16 \
  --fastq /path/to/reads.fastq \
  -m TTAGGG
```

**Assembly-only comparison** (benchmarking only):

```bash
python3 -m taco -g 12m -t 16 \
  --fastq /path/to/reads.fastq \
  -m TTAGGG \
  --assembly-only
```

**With Merqury QV scoring:**

```bash
python3 -m taco -g 12m -t 16 \
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
| `--fasta` | External pre-assembled FASTA to include in comparison |
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
├── run_taco                # Shell wrapper (runs python3 -m taco)
├── taco/                   # Python package
│   ├── __init__.py         # Package metadata (v1.0.0)
│   ├── __main__.py         # Entry point (python3 -m taco)
│   ├── cli.py              # Argument parsing
│   ├── pipeline.py         # Pipeline runner, logging, benchmarking
│   ├── steps.py            # All 18 step implementations
│   ├── utils.py            # Shared utilities and FASTA I/O
│   ├── telomere_detect.py  # Hybrid telomere detection engine
│   ├── telomere_pool.py    # Telomere pool classification
│   ├── clustering.py       # Minimap2-based contig clustering
│   ├── backbone.py         # Backbone selection and scoring
│   └── reporting.py        # Final report generation
├── taco-env.yml            # Conda environment specification
├── INSTALLATION.md         # Installation guide
├── README.md
├── LICENSE                 # MIT License
└── .gitignore
```

## Troubleshooting

**Canu reports `master +XX changes`:** You are using a development build. Install a stable Canu release before benchmarking.

**IPA or Peregrine skipped:** These assemblers only support certain platforms. IPA requires PacBio HiFi; Peregrine does not support Nanopore. Use `--platform` to match your data type.

**Telomere motif appears incorrect:** Do not assume the same motif for all fungi. Choose a motif appropriate for the target organism. TACO's built-in motif families cover canonical TTAGGG, budding yeast TG1-3, and Candida repeats.

**`TACO.sh: command not found`:** Add the TACO directory to your PATH or run with the full path.

**Missing Python modules:** TACO uses only the Python standard library. If you see import errors, ensure Python >= 3.8 is installed and the `taco/` directory is alongside `TACO.sh`.

**Merqury not working:** Merqury is optional. Install with `conda install -c bioconda merqury meryl` or use `--no-merqury` to skip.

## Citation

If you use TACO in a publication, please cite the software and archive the exact release used for reproducibility (e.g., via Zenodo).

TACO was developed at the Grainger Bioinformatics Center, Field Museum of Natural History, Chicago, Illinois, USA.

## License

TACO is released under the [MIT License](LICENSE).
