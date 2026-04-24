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

**Assemblers (via conda):** Canu, NextDenovo, Flye, Hifiasm, Peregrine, IPA, LJA, Raven

**HiCanu:** The conda environment includes `canu` and `openjdk>=11` to provide a working Java runtime and avoid the `undefined symbol: JLI_StringDup` error. If you still get Java errors (e.g. from a bioconda dev build), download a stable binary from https://github.com/marbl/canu/releases and place it on PATH. If canu is missing or fails, Step 1 is skipped and all other assemblers continue normally.

**Analysis tools (via conda):** BUSCO, QUAST, Minimap2, Seqtk, BWA, Samtools, purge_dups

**Polishing tools (via conda):** NextPolish2, yak, Racon, Medaka

**QV scoring (via conda):** Merqury, Meryl (auto-enabled for HiFi data; builds reads.meryl from input reads if no pre-built database exists)

**Optional manual install:** MBG (Multiplex de Bruijn Graph, HiFi-only). Build from source: `git clone https://github.com/maickrau/MBG && cd MBG && make`. Place the `MBG` binary on PATH.

## Sequencing Platform Support

TACO supports three sequencing platforms via the `--platform` flag:

| Platform | Flag | Assemblers Used |
|----------|------|-----------------|
| PacBio HiFi | `--platform pacbio-hifi` (default) | canu, nextDenovo, flye, hifiasm, peregrine, ipa, lja, mbg, raven |
| Oxford Nanopore | `--platform nanopore` | canu, nextDenovo, flye, hifiasm (UL mode), raven |
| PacBio CLR | `--platform pacbio` | canu, nextDenovo, flye, peregrine, raven |

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

```bash
mkdir -p my_project && cd my_project

# Fungal genome (HiFi reads)
taco -g 12m -t 16 --fastq /path/to/reads.fastq.gz --taxon fungal

# Plant genome (HiFi reads)
taco -g 500m -t 32 --fastq /path/to/reads.fastq.gz --taxon plant

# Vertebrate genome (ONT reads)
taco -g 2.5g -t 32 --fastq /path/to/reads.fastq.gz --taxon vertebrate --platform nanopore

# Assembly-only comparison (no refinement)
taco -g 12m -t 16 --fastq /path/to/reads.fastq.gz --taxon fungal --assembly-only
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
For HiFi data, Merqury is auto-enabled if `merqury.sh` and `meryl` are installed. The reads.meryl database is built automatically. If you have a pre-built database, provide it with `--merqury-db path/to/reads.meryl`. Disable with `--no-merqury`.

**NextPolish2 requires samtools:**
NextPolish2 v0.2+ needs a sorted BAM file. TACO maps reads with minimap2 and sorts with samtools. Ensure `samtools` is installed in the conda env.

**Canu reports `master +XX changes` or Step 1 fails with a Java error:**
The conda environment includes `openjdk>=11`. If you still see Java errors, download a stable canu binary from https://github.com/marbl/canu/releases. TACO detects dev builds and warns you. If canu fails, the pipeline continues with the remaining assemblers.
