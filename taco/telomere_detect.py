"""Hybrid telomere detection engine for TACO."""
import os
import sys
import re
import csv
import argparse
import subprocess
from collections import defaultdict
import math

try:
    from taco.utils import read_fasta, write_fasta, revcomp
except ImportError:
    try:
        from utils import read_fasta, write_fasta, revcomp
    except ImportError:
        pass


MOTIF_FAMILIES = {
    "canonical": {
        "motif": "TTAGGG",
        "revcomp": "CCCTAA",
        "description": "Vertebrate/filamentous fungal canonical repeat",
        "regex_fwd": r"(TTAGGG){2,}",
        "regex_rev": r"(CCCTAA){2,}",
    },
    "budding_yeast": {
        "motif": "TG1-3",
        "description": "Budding yeast degenerate TG(1-3)/C(1-3)A",
        "regex_fwd": r"(TG{1,3}){3,}",
        "regex_rev": r"(C{1,3}A){3,}",
    },
    "candida": {
        "motif": "ACGGATGTCTAACTTCTTGGTGT",
        "description": "Candida 23-bp telomere repeat",
        "regex_fwd": r"(ACGGATGTCTAACTTCTTGGTGT){2,}",
        "regex_rev": r"(ACACCAAGAAGTTTAGACATCCGT){2,}",
    },
}

SCORING_WEIGHTS = {
    "density": 0.40,
    "longest_run": 0.30,
    "distance_to_end": 0.20,
    "covered_bp": 0.10,
}

THRESHOLDS = {
    "strict_t2t": 0.25,
    "single_tel_strong": 0.25,
    "telomere_supported": 0.08,
}


def extract_ends(seq, window=5000):
    """Extract terminal sequences from contig.

    Args:
        seq: DNA sequence string
        window: Length of terminal window to extract (default 5000)

    Returns:
        tuple: (left_end, right_end)
    """
    left_end = seq[:window] if len(seq) >= window else seq
    right_end = seq[-window:] if len(seq) >= window else seq
    return left_end, right_end


def build_regex_for_motif(motif):
    """Build compiled regex for forward and reverse complement.

    Args:
        motif: Motif string (can be degenerate like TG{1,3})

    Returns:
        tuple: (compiled_fwd_regex, compiled_rev_regex)
    """
    try:
        fwd_regex = re.compile(motif, re.IGNORECASE)
    except:
        fwd_regex = None

    try:
        rev_seq = revcomp(motif)
        rev_regex = re.compile(rev_seq, re.IGNORECASE)
    except:
        rev_regex = None

    return fwd_regex, rev_regex


def discover_motifs_python(seqs_dict, end_window=5000, kmer_min=4, kmer_max=30, top_n=10):
    """Discover telomeric motifs by k-mer frequency analysis of contig ends.

    Args:
        seqs_dict: dict of {name: sequence}
        end_window: Length of contig ends to analyze (default 5000)
        kmer_min: Minimum k-mer size (default 4)
        kmer_max: Maximum k-mer size (default 30)
        top_n: Number of top motifs to return (default 10)

    Returns:
        list: Discovered motif strings
    """
    kmer_counts = defaultdict(int)

    for name, seq in seqs_dict.items():
        # Extract ends
        left_end = seq[:end_window]
        right_end = seq[-end_window:] if len(seq) > end_window else seq

        # Count k-mers
        for kmer_size in range(kmer_min, kmer_max + 1):
            for end in [left_end, right_end]:
                for i in range(len(end) - kmer_size + 1):
                    kmer = end[i:i+kmer_size].upper()
                    kmer_counts[kmer] += 1

    # Sort by frequency and return top n
    sorted_kmers = sorted(kmer_counts.items(), key=lambda x: x[1], reverse=True)
    return [kmer for kmer, count in sorted_kmers[:top_n]]


def discover_motifs_tidk(fasta_path, threads=1):
    """Discover motifs using tidk if available.

    Args:
        fasta_path: Path to FASTA file
        threads: Number of threads (default 1)

    Returns:
        list: Motif strings discovered by tidk, or empty list if tidk unavailable
    """
    try:
        result = subprocess.run(
            ["tidk", "search", "--fasta", fasta_path, "--threads", str(threads)],
            capture_output=True,
            text=True,
            timeout=300
        )
        motifs = []
        for line in result.stdout.split('\n'):
            if line.strip() and not line.startswith('#'):
                parts = line.split()
                if parts:
                    motifs.append(parts[0])
        return motifs
    except:
        return []


def score_end(end_seq, motifs_regexes, score_window=500):
    """Score a sequence end for telomeric content.

    Args:
        end_seq: Terminal sequence to score
        motifs_regexes: List of (fwd_regex, rev_regex) tuples
        score_window: Window size for scoring (default 500)

    Returns:
        dict: Keys are density, longest_run, distance_to_end, covered_bp, raw_score
    """
    window = end_seq[:score_window] if len(end_seq) >= score_window else end_seq

    best_matches = []
    total_bp_covered = 0
    longest_run = 0

    for fwd_regex, rev_regex in motifs_regexes:
        if fwd_regex is None:
            continue

        for match in fwd_regex.finditer(window):
            best_matches.append((match.start(), match.end()))
            total_bp_covered += match.end() - match.start()
            run_length = match.end() - match.start()
            if run_length > longest_run:
                longest_run = run_length

        if rev_regex is not None:
            for match in rev_regex.finditer(window):
                best_matches.append((match.start(), match.end()))
                total_bp_covered += match.end() - match.start()
                run_length = match.end() - match.start()
                if run_length > longest_run:
                    longest_run = run_length

    # Merge overlapping matches
    if best_matches:
        best_matches.sort()
        merged = [best_matches[0]]
        for start, end in best_matches[1:]:
            if start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        total_bp_covered = sum(e - s for s, e in merged)

    window_len = len(window)
    if window_len == 0:
        return {
            "density": 0.0,
            "longest_run": 0.0,
            "distance_to_end": 0.0,
            "covered_bp": 0.0,
            "raw_score": 0.0,
        }

    density = total_bp_covered / window_len
    longest_run_norm = longest_run / window_len
    distance_to_end = 1.0 - (0 if not best_matches else best_matches[-1][1]) / window_len
    distance_to_end = max(0, min(1.0, distance_to_end))
    covered_bp_norm = total_bp_covered / window_len

    raw_score = (
        SCORING_WEIGHTS["density"] * density +
        SCORING_WEIGHTS["longest_run"] * longest_run_norm +
        SCORING_WEIGHTS["distance_to_end"] * distance_to_end +
        SCORING_WEIGHTS["covered_bp"] * covered_bp_norm
    )

    return {
        "density": density,
        "longest_run": longest_run_norm,
        "distance_to_end": distance_to_end,
        "covered_bp": covered_bp_norm,
        "raw_score": raw_score,
    }


def classify_contig(left_score, right_score):
    """Classify contig based on telomere scores.

    Args:
        left_score: Score of left telomere (0-1)
        right_score: Score of right telomere (0-1)

    Returns:
        str: Classification category
    """
    if left_score >= THRESHOLDS["strict_t2t"] and right_score >= THRESHOLDS["strict_t2t"]:
        return "strict_t2t"
    elif max(left_score, right_score) >= THRESHOLDS["single_tel_strong"]:
        return "single_tel_strong"
    elif max(left_score, right_score) >= THRESHOLDS["telomere_supported"]:
        return "telomere_supported"
    else:
        return "none"


def detect_telomeres(fasta_path, mode="hybrid", user_motif=None, end_window=5000,
                     score_window=500, kmer_min=4, kmer_max=30, threads=1):
    """Main telomere detection function.

    Args:
        fasta_path: Path to FASTA file
        mode: Detection mode - 'known', 'auto', 'hybrid'
        user_motif: User-provided motif string (optional)
        end_window: Length of terminal windows to analyze
        score_window: Window size for scoring
        kmer_min: Minimum k-mer size for discovery
        kmer_max: Maximum k-mer size for discovery
        threads: Number of threads for tidk

    Returns:
        list: Dicts with contig, left_score, right_score, classification
    """
    seqs = read_fasta(fasta_path)

    # Gather motifs based on mode
    motifs = []
    if mode == "known" and user_motif:
        motifs = [user_motif]
    elif mode == "auto":
        motifs = discover_motifs_python(seqs, end_window, kmer_min, kmer_max)
    elif mode == "hybrid":
        if user_motif:
            motifs = [user_motif]
        discovered = discover_motifs_python(seqs, end_window, kmer_min, kmer_max, top_n=3)
        motifs.extend(discovered)
        for fam_name, fam_data in MOTIF_FAMILIES.items():
            motifs.append(fam_data["regex_fwd"])

    # Build regexes
    motifs_regexes = []
    for motif in motifs:
        if isinstance(motif, str):
            fwd_regex, rev_regex = build_regex_for_motif(motif)
            if fwd_regex:
                motifs_regexes.append((fwd_regex, rev_regex))

    # Score each contig
    results = []
    for contig_name, seq in seqs.items():
        left_end, right_end = extract_ends(seq, end_window)

        left_result = score_end(left_end, motifs_regexes, score_window)
        right_result = score_end(right_end, motifs_regexes, score_window)

        left_score = left_result["raw_score"]
        right_score = right_result["raw_score"]
        classification = classify_contig(left_score, right_score)

        results.append({
            "contig": contig_name,
            "left_score": left_score,
            "right_score": right_score,
            "classification": classification,
        })

    return results


def write_detection_outputs(results, fasta_path, out_prefix):
    """Write telomere detection outputs.

    Args:
        results: List of detection result dicts
        fasta_path: Path to input FASTA
        out_prefix: Prefix for output files
    """
    seqs = read_fasta(fasta_path)

    # Write .telomere_end_scores.tsv
    scores_path = f"{out_prefix}.telomere_end_scores.tsv"
    with open(scores_path, 'w') as f:
        f.write("contig\tleft_score\tright_score\tclassification\n")
        for result in results:
            f.write(f"{result['contig']}\t{result['left_score']:.4f}\t{result['right_score']:.4f}\t{result['classification']}\n")

    # Write .telo_metrics.tsv
    metrics_path = f"{out_prefix}.telo_metrics.tsv"
    counts = defaultdict(int)
    for result in results:
        counts[result["classification"]] += 1

    with open(metrics_path, 'w') as f:
        f.write("metric\tvalue\n")
        f.write(f"strict_t2t\t{counts['strict_t2t']}\n")
        f.write(f"single_tel_strong\t{counts['single_tel_strong']}\n")
        f.write(f"telomere_supported\t{counts['telomere_supported']}\n")
        f.write(f"total_contigs\t{len(results)}\n")

    # Write .telo.fasta and .telo.list
    telo_seqs = {}
    telo_list = []
    for result in results:
        if result["classification"] != "none":
            contig_name = result["contig"]
            telo_seqs[contig_name] = seqs[contig_name]
            telo_list.append(contig_name)

    telo_fasta_path = f"{out_prefix}.telo.fasta"
    write_fasta(telo_seqs, telo_fasta_path)

    telo_list_path = f"{out_prefix}.telo.list"
    with open(telo_list_path, 'w') as f:
        for name in telo_list:
            f.write(f"{name}\n")


def main():
    """Command-line interface for telomere detection."""
    parser = argparse.ArgumentParser(description="Hybrid telomere detection")
    parser.add_argument("--fasta", required=True, help="Input FASTA file")
    parser.add_argument("--mode", default="hybrid", choices=["known", "auto", "hybrid"],
                        help="Detection mode")
    parser.add_argument("--out-prefix", required=True, help="Output file prefix")
    parser.add_argument("--motif", help="User-provided motif")
    parser.add_argument("--end-window", type=int, default=5000, help="Terminal window size")
    parser.add_argument("--score-window", type=int, default=500, help="Scoring window size")
    parser.add_argument("--kmer-min", type=int, default=4, help="Minimum k-mer size")
    parser.add_argument("--kmer-max", type=int, default=30, help="Maximum k-mer size")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads")

    args = parser.parse_args()

    results = detect_telomeres(
        args.fasta,
        mode=args.mode,
        user_motif=args.motif,
        end_window=args.end_window,
        score_window=args.score_window,
        kmer_min=args.kmer_min,
        kmer_max=args.kmer_max,
        threads=args.threads,
    )

    write_detection_outputs(results, args.fasta, args.out_prefix)

    print(f"Telomere detection complete. Output prefix: {args.out_prefix}")


if __name__ == "__main__":
    main()
