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

**Assemblers (via conda):** NextDenovo, Flye, Hifiasm

**HiCanu:** The conda environment includes `canu` and `openjdk>=11` to provide a working Java runtime and avoid the `undefined symbol: JLI_StringDup` error. If you still get Java errors (e.g. from a bioconda dev build), download a stable binary from https://github.com/marbl/canu/releases and place it on PATH. If canu is missing or fails, Step 1 is skipped and all other assemblers continue normally.

**Analysis tools (via conda):** BUSCO, QUAST, Minimap2, Funannotate, Seqtk, BWA, Samtools

**Optional tools (via conda):** Merqury + Meryl (assembly QV scoring), tidk (telomere repeat discovery)

## Manual Installation of Peregrine and IPA

Peregrine and IPA are not reliably available through current bioconda channels and may require manual installation.

### Peregrine

```bash
# Clone and build from source
git clone https://github.com/cschin/peregrine-2021.git
cd peregrine-2021
# Follow the build instructions in the repository
```

Ensure `pg_asm` is on your PATH after installation. If Peregrine is not installed, TACO will skip Step 3 with a warning.

### IPA (PacBio Improved Phased Assembler)

```bash
# Try conda first
conda install -c bioconda pbipa

# Or install from source
# https://github.com/PacificBiosciences/pbipa
```

IPA only supports PacBio HiFi reads. If not installed, TACO will skip Step 4 with a warning.

## Running TACO

After installation, simply use the `taco` command:

```bash
mkdir -p my_project && cd my_project
taco -g 12m -t 16 --fastq /path/to/reads.fastq -m TTAGGG
```

**Alternative (without pip install):** Use the shell wrapper which sets PYTHONPATH automatically:

```bash
~/opt/TACO/run_taco -g 12m -t 16 --fastq /path/to/reads.fastq
```

## Sequencing Platform Support

TACO supports three sequencing platforms via the `--platform` flag:

| Platform | Flag | Assemblers Used |
|----------|------|-----------------|
| PacBio HiFi | `--platform pacbio-hifi` (default) | All 6 assemblers |
| Oxford Nanopore | `--platform nanopore` | canu, nextDenovo, flye, hifiasm |
| PacBio CLR | `--platform pacbio` | canu, nextDenovo, peregrine, flye, hifiasm |

Incompatible or missing assemblers are automatically skipped with a warning.

## Troubleshooting

**`taco: command not found` after pip install:**
Make sure the taco conda environment is activated: `conda activate taco`

**`No module named taco`:**
Run `pip install -e .` from the TACO repository directory, or use `./run_taco` instead.

**Missing Python modules:**
TACO uses only the Python standard library. If you see import errors, ensure Python >= 3.8 is installed and run `pip install -e .` again.

**Merqury not working:**
Merqury is optional. Install with `conda install -c bioconda merqury meryl` or use `--no-merqury` to skip.

**Canu reports `master +XX changes` or Step 1 fails with a Java error:**
The conda environment now includes `openjdk>=11` to provide a working Java runtime. If you still see this error, the bioconda canu package may be a dev build. Download a stable binary from https://github.com/marbl/canu/releases and place it on PATH. TACO detects dev builds and warns you. If canu fails, the pipeline continues with the remaining assemblers.
