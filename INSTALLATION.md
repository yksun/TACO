# TACO Installation Guide

## Prerequisites

TACO requires a Unix-like operating system (Linux or macOS) with:

- Bash >= 4.0
- Python >= 3.8
- Conda or Micromamba

## Quick Install

```bash
# Clone the repository
git clone https://github.com/ysun-fieldmuseum/TACO.git
cd TACO

# Create the conda environment
conda env create -f taco-env.yml

# Activate the environment
conda activate taco

# Verify installation
python3 -m taco --help
```

## Conda Environment

TACO bundles all dependencies in a single conda environment defined in
`taco-env.yml`. The environment includes:

**Assemblers:** HiCanu (via canu), NextDenovo, Peregrine, IPA (PacBio), Flye, Hifiasm

**Analysis tools:** BUSCO, QUAST, Minimap2, Funannotate, Seqtk, BWA, Samtools

**Optional tools:** Merqury + Meryl (assembly QV scoring), tidk (improves telomere auto-detection)

## Manual Installation

If you prefer not to use the bundled environment:

1. Install Python >= 3.8 (standard library only; no pip packages needed)
2. Install each assembler following its own documentation
3. Install analysis tools: busco, quast, minimap2, funannotate, seqtk, bwa, samtools
4. Ensure all tools are available on your PATH
5. Run `python3 -m taco --help` to verify

## Sequencing Platform Support

TACO supports three sequencing platforms via the `--platform` flag:

| Platform | Flag | Assemblers Used |
|----------|------|-----------------|
| PacBio HiFi | `--platform pacbio-hifi` (default) | All 6 assemblers |
| Oxford Nanopore | `--platform nanopore` | canu, nextDenovo, flye, hifiasm |
| PacBio CLR | `--platform pacbio` | canu, nextDenovo, peregrine, flye, hifiasm |

IPA only supports PacBio HiFi reads. Peregrine does not support Nanopore reads.
Incompatible assemblers are automatically skipped with a warning.

## Running TACO

TACO can be run in two ways:

```bash
# Option 1: Python module (recommended)
python3 -m taco -g 12m -t 16 --fastq reads.fastq

# Option 2: Shell wrapper
./run_taco -g 12m -t 16 --fastq reads.fastq
```

To run from any directory, add the TACO directory to PYTHONPATH:

```bash
echo 'export PYTHONPATH="/path/to/TACO:$PYTHONPATH"' >> ~/.bashrc
source ~/.bashrc
```

## Project Structure

```
TACO/
├── run_taco                # Shell wrapper (runs python3 -m taco)
├── taco/                   # Python package
│   ├── __init__.py         # Package metadata (v1.0.0)
│   ├── __main__.py         # Entry point (python3 -m taco)
│   ├── cli.py              # Argument parsing
│   ├── pipeline.py         # Pipeline runner
│   ├── steps.py            # All 18 step implementations
│   ├── utils.py            # Shared utilities and FASTA I/O
│   ├── telomere_detect.py  # Hybrid telomere detection engine
│   ├── telomere_pool.py    # Telomere pool classification
│   ├── clustering.py       # Minimap2-based contig clustering
│   ├── backbone.py         # Backbone selection and scoring
│   └── reporting.py        # Final report generation
├── taco-env.yml            # Conda environment specification
├── INSTALLATION.md
├── README.md
├── LICENSE
└── .gitignore
```

## Troubleshooting

**"Command not found" errors:**
Ensure the taco conda environment is activated: `conda activate taco`

**Permission denied on run_taco:**
Make the wrapper executable: `chmod +x run_taco`, or use `python3 -m taco` directly.

**Missing Python modules:**
TACO uses only the Python standard library. If you see import errors, ensure
Python >= 3.8 is installed and the `taco/` directory is alongside `TACO.sh`.

**Merqury not working:**
Merqury is optional. If not installed, TACO will skip Merqury metrics.
To install: `conda install -c bioconda merqury meryl`

**Canu reports `master +XX changes`:**
You are using a development build. Install a stable Canu release (e.g., v2.2 or v2.3)
before running benchmarks.
