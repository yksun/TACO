"""Utility functions for TACO bioinformatics pipeline."""
import os
import sys
import csv
import re
from collections import defaultdict


# Canonical assembler list — all downstream code should import this instead
# of hardcoding assembler names.  "reference" is included for optional
# user-provided reference assembly comparison.
ALL_ASSEMBLERS = [
    "canu", "reference", "flye", "ipa", "nextDenovo",
    "peregrine", "hifiasm", "lja", "mbg", "raven",
]

ASSEMBLER_PLATFORMS = {
    "canu": {"pacbio-hifi": "-pacbio-hifi", "nanopore": "-nanopore", "pacbio": "-pacbio"},
    "nextDenovo": {"pacbio-hifi": "hifi", "nanopore": "ont", "pacbio": "clr"},
    "peregrine": {"pacbio-hifi": True, "nanopore": False, "pacbio": True},
    "ipa": {"pacbio-hifi": True, "nanopore": False, "pacbio": False},
    "flye": {"pacbio-hifi": "--pacbio-hifi", "nanopore": "--nano-hq", "pacbio": "--pacbio-raw"},
    "hifiasm": {"pacbio-hifi": True, "nanopore": False, "pacbio": False},
    # hifiasm: HiFi-only as primary input. ONT ultra-long reads can supplement
    # HiFi via --ul flag but require HiFi as primary — use --enable-hifiasm-ul
    # with paired HiFi+ONT data (not implemented in TACO auto-mode yet).
    "lja": {"pacbio-hifi": True, "nanopore": False, "pacbio": False},
    # LJA (La Jolla Assembler): HiFi-only, produces very contiguous assemblies
    "mbg": {"pacbio-hifi": True, "nanopore": False, "pacbio": False},
    # MBG (Multiplex de Bruijn Graph): HiFi-only, minimizer-based
    "raven": {"pacbio-hifi": True, "nanopore": True, "pacbio": True},
    # raven: supports all long-read platforms, fast OLC assembler
}


def get_assembler_flag(assembler, platform):
    """Get the assembler-specific flag for a platform."""
    return ASSEMBLER_PLATFORMS.get(assembler, {}).get(platform)


def is_assembler_compatible(assembler, platform):
    """Check if an assembler supports a platform."""
    flag = get_assembler_flag(assembler, platform)
    return flag is not None and flag is not False


def read_fasta(path):
    """Parse FASTA file and return dict of {name: sequence}.

    Args:
        path: Path to FASTA file

    Returns:
        dict: {header (first word after >): concatenated sequence}
    """
    seqs = {}
    current_name = None
    current_seq = []

    if not os.path.exists(path):
        raise FileNotFoundError(f"FASTA file not found: {path}")

    with open(path, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_name is not None:
                    seqs[current_name] = ''.join(current_seq)
                # Extract header (first word after >)
                header = line[1:].split()[0]
                current_name = header
                current_seq = []
            else:
                current_seq.append(line)

        # Save last sequence
        if current_name is not None:
            seqs[current_name] = ''.join(current_seq)

    return seqs


def write_fasta(seqs_dict, path, wrap=60):
    """Write dict to FASTA with line wrapping.

    Args:
        seqs_dict: dict of {name: sequence}
        path: Output FASTA path
        wrap: Line wrap length (default 60)
    """
    with open(path, 'w') as f:
        for name, seq in seqs_dict.items():
            f.write(f">{name}\n")
            for i in range(0, len(seq), wrap):
                f.write(seq[i:i+wrap] + '\n')


def count_fasta(path):
    """Count sequences in FASTA file by counting > lines.

    Args:
        path: Path to FASTA file

    Returns:
        int: Number of sequences
    """
    count = 0
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count


def filter_fasta_by_ids(input_path, ids_set, output_path):
    """Write subset of FASTA matching ids.

    Args:
        input_path: Input FASTA path
        ids_set: set of sequence IDs to keep
        output_path: Output FASTA path
    """
    seqs = read_fasta(input_path)
    filtered = {name: seq for name, seq in seqs.items() if name in ids_set}
    write_fasta(filtered, output_path)


def rename_and_sort_fasta(input_path, output_path, prefix):
    """Sort contigs by length descending, rename as {prefix}_1, {prefix}_2, etc., wrap at 60bp.

    Args:
        input_path: Input FASTA path
        output_path: Output FASTA path
        prefix: Prefix for new contig names
    """
    seqs = read_fasta(input_path)
    # Sort by length descending
    sorted_seqs = sorted(seqs.items(), key=lambda x: len(x[1]), reverse=True)

    renamed = {}
    for idx, (orig_name, seq) in enumerate(sorted_seqs, 1):
        renamed[f"{prefix}_{idx}"] = seq

    write_fasta(renamed, output_path, wrap=60)


def revcomp(seq):
    """Return reverse complement of DNA sequence.

    Args:
        seq: DNA sequence string

    Returns:
        str: Reverse complement
    """
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                  'N': 'N', 'n': 'n'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


def read_tsv_column(path, col_name):
    """Read a specific column from a TSV file.

    Args:
        path: Path to TSV file
        col_name: Column name to extract

    Returns:
        list: Values from the column
    """
    values = []
    with open(path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if reader.fieldnames is None or col_name not in reader.fieldnames:
            raise ValueError(f"Column '{col_name}' not found in {path}")
        for row in reader:
            values.append(row[col_name])
    return values


def merge_csv_files(file_list, output_path):
    """Merge multiple CSVs with same header format.

    Args:
        file_list: list of CSV file paths
        output_path: Output CSV path
    """
    all_rows = []
    fieldnames = None

    for file_path in file_list:
        if not os.path.exists(file_path):
            continue
        with open(file_path, 'r') as f:
            reader = csv.DictReader(f)
            if fieldnames is None:
                fieldnames = reader.fieldnames
            for row in reader:
                all_rows.append(row)

    if fieldnames is None:
        fieldnames = []

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_rows)
