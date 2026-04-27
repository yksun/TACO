# TACO Installation Guide

## Prerequisites

TACO requires a Unix-like operating system (Linux or macOS) with:

- Python >= 3.8
- Conda or Micromamba

## Quick Install

```bash
# Clone the repository
git clone https://github.com/yksun/TACO.git
cd TACO

# Create the conda environment (includes all analysis tools)
conda env create -f taco-env.yml
conda activate taco

# Install TACO as a command
pip install -e .

# Verify
taco --help
```

After installation, the `taco` command is available anywhere within the conda environment.

## What Gets Installed

The conda environment provides all external bioinformatics tools. The `pip install -e .` step registers the `taco` command so you can run it from any directory.

**Assemblers (via conda):** Canu, NextDenovo, Flye, Hifiasm, Peregrine, IPA, LJA, Raven (`raven-assembler`)

**HiCanu:** The conda environment includes `canu` and `openjdk>=11` to provide a working Java runtime and avoid the `undefined symbol: JLI_StringDup` error. If you still get Java errors (e.g. from a bioconda dev build), download a stable binary from https://github.com/marbl/canu/releases and place it on PATH. If canu is missing or fails, Step 1 is skipped and all other assemblers continue normally.

**Analysis tools (via conda):** BUSCO, QUAST, Minimap2, Seqtk, BWA, Samtools, purge_dups

**Polishing tools (via conda):** NextPolish2, yak, Racon, Medaka

**QV scoring (via conda):** Merqury, Meryl (auto-enabled for all platforms when installed; builds a reads `.meryl` database from input reads with `--merqury-k auto` by default. QV is most accurate with PacBio HiFi or Illumina; ONT/CLR QV may underestimate true quality, and relative ranking/completeness should be interpreted cautiously)

**Additional HiFi assembler:** MBG (Multiplex de Bruijn Graph, HiFi-only) is included in `taco-env.yml` via Bioconda. For an existing environment, install with `mamba install -c conda-forge -c bioconda mbg`. TACO passes MBG's required odd `-k` k-mer size with a default of `1501`; set `TACO_MBG_K` to override. If the package is unavailable on your platform, build from source: `git clone https://github.com/maickrau/MBG && cd MBG && make`, then place the `MBG` binary on PATH.

## Sequencing Platform Support

TACO supports three sequencing platforms via the `--platform` flag:

| Platform | Flag | Assemblers Used |
|----------|------|-----------------|
| PacBio HiFi | `--platform pacbio-hifi` (default) | canu, nextDenovo, flye, hifiasm, peregrine, ipa, lja, mbg, raven (up to 9) |
| Oxford Nanopore | `--platform nanopore` | canu, nextDenovo, flye, raven (4 assemblers) |
| PacBio CLR | `--platform pacbio` | canu, nextDenovo, flye, peregrine, raven (5 assemblers) |

Incompatible or missing assemblers are automatically skipped with a warning.

## Taxon-Aware Defaults

TACO automatically selects appropriate BUSCO lineage, scoring weights, and telomere motifs based on `--taxon`:

| Taxon | Default BUSCO | Telomere Motifs |
|-------|--------------|-----------------|
| `--taxon fungal` | ascomycota_odb10 | TTAGGG + TG1-3 + Candida |
| `--taxon plant` | embryophyta_odb10 | TTTAGGG |
| `--taxon vertebrate` | vertebrata_odb10 | TTAGGG |
| `--taxon insect` | insecta_odb10 | TTAGG |
| `--taxon other` | requires `--busco` | all known motif families |

Override the BUSCO lineage with `--busco <lineage_name>` if your organism needs a more specific database.

## Running TACO

TACO runs 17 public steps (0-16). Step 0 (Input QC) always runs first — it validates the FASTQ, estimates coverage, and logs compatible assemblers. In full mode (Steps 0-15), TACO assembles, runs combined assembly QC/comparison, builds the telomere pool, selects a backbone, refines it, and produces a final assembly. In assembly-only mode (Steps 0-11 and 16), TACO benchmarks all assemblers without refinement.

```bash
mkdir -p my_project && cd my_project

# Full refinement: fungal genome (HiFi reads)
taco -g 12m -t 16 --fastq /path/to/reads.fastq.gz --taxon fungal

# Full refinement: plant genome (HiFi reads)
taco -g 500m -t 32 --fastq /path/to/reads.fastq.gz --taxon plant

# Full refinement: vertebrate genome (ONT reads)
taco -g 2.5g -t 32 --fastq /path/to/reads.fastq.gz --taxon vertebrate --platform nanopore

# Assembly-only comparison (benchmarking only, no refinement)
taco -g 12m -t 16 --fastq /path/to/reads.fastq.gz --taxon fungal --assembly-only

# Optional step timing/provenance logs
taco -g 12m -t 16 --fastq /path/to/reads.fastq.gz --taxon fungal --benchmark

# Optional publication-grade input checksums with benchmark logs
TACO_BENCHMARK_SHA256=1 taco -g 12m -t 16 --fastq /path/to/reads.fastq.gz --taxon fungal --benchmark

# Resume from a specific step (e.g., rerun from refinement)
taco -g 12m -t 16 --fastq /path/to/reads.fastq.gz --taxon fungal -s 13-15
```

**Alternative (without pip install):** Use the shell wrapper which sets PYTHONPATH automatically:

```bash
~/opt/TACO/run_taco -g 12m -t 16 --fastq /path/to/reads.fastq.gz --taxon fungal
```

## Troubleshooting

**`taco: command not found` after pip install:**
Make sure the taco conda environment is activated: `conda activate taco`

**`No module named taco`:**
Run `pip install -e .` from the TACO repository directory, or use `./run_taco` instead.

**Missing Python modules:**
TACO uses only the Python standard library. If you see import errors, ensure Python >= 3.8 is installed and run `pip install -e .` again.

**Merqury not working:**
Merqury is auto-enabled for all platforms when `merqury.sh` and `meryl` are installed. The reads `.meryl` database is built automatically from input reads. `--merqury-k auto` is the default and uses Merqury `best_k.sh` when available, with a genome-size fallback. For Nanopore/PacBio CLR data, QV values may underestimate true quality (TACO logs a warning), and Merqury completeness/relative ranking should be interpreted cautiously. Provide a pre-built database with `--merqury-db path/to/reads.meryl`. Disable with `--no-merqury`.

**NextPolish2 requires samtools:**
NextPolish2 v0.2+ needs a sorted BAM file. TACO maps reads with minimap2 and sorts with samtools. Ensure `samtools` is installed in the conda env.

**Canu reports `master +XX changes` or Step 1 fails with a Java error:**
The conda environment includes `openjdk>=11`. If you still see Java errors, download a stable canu binary from https://github.com/marbl/canu/releases. TACO detects dev builds and warns you. If canu fails, the pipeline continues with the remaining assemblers.
