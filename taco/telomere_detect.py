"""Hybrid telomere detection engine for TACO.

Faithfully implements the same scoring algorithm as TACO.sh:
  - De novo k-mer discovery with enrichment scoring and canonical rotation
  - MOTIF_FAMILIES (canonical TTAGGG, budding yeast TG1-3, Candida 23-bp)
  - Per-end composite scoring: density*0.40 + longest_run*0.30 +
    distance_to_end*0.20 + covered_bp*0.10
  - Three-tier classification: strict_t2t, single_tel_strong, telomere_supported
"""
import os
import sys
import re
import csv
import argparse
import subprocess
import shutil
import tempfile
from collections import Counter, defaultdict
from pathlib import Path


# ── Motif family definitions (matches TACO.sh exactly) ──────────────────────
MOTIF_FAMILIES = {
    "canonical": {
        "forward": ["TTAGGG"],
        "regex_fwd": r"(TTAGGG){2,}",
        "regex_rev": r"(CCCTAA){2,}",
    },
    "budding_yeast": {
        "forward": [],
        "regex_fwd": r"(TG{1,3}){3,}",
        "regex_rev": r"(C{1,3}A){3,}",
    },
    "candida": {
        "forward": ["ACGGATGTCTAACTTCTTGGTGT"],
        "regex_fwd": r"(ACGGATGTCTAACTTCTTGGTGT){2,}",
        "regex_rev": r"(ACACCAAGAAGTTAGACATCCGT){2,}",
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


# ── Taxon-aware telomere presets ────────────────────────────────────────────
# Each preset defines which motif families to prioritize and default behavior.
# --motif overrides these priors when the user knows the exact repeat.
TAXON_PRESETS = {
    "vertebrate": {
        "primary_families": ["canonical"],
        "description": "TTAGGG — highly conserved in vertebrates",
    },
    "animal": {
        "primary_families": ["canonical"],
        "description": "TTAGGG — strongest prior for vertebrates, less certain for distant metazoans",
    },
    "plant": {
        "primary_families": ["plant"],
        "description": "TTTAGGG — common plant telomere repeat",
    },
    "insect": {
        "primary_families": ["insect"],
        "description": "TTAGG — common insect telomere repeat",
    },
    "fungal": {
        "primary_families": ["canonical", "budding_yeast", "candida"],
        "description": "Diverse fungal telomeres — uses all built-in families",
    },
    "other": {
        "primary_families": ["canonical", "budding_yeast", "candida", "plant", "insect"],
        "description": "Unknown taxon — uses all motif families plus de novo discovery",
    },
}


# Extended MOTIF_FAMILIES with plant and insect entries
MOTIF_FAMILIES["plant"] = {
    "forward": ["TTTAGGG"],
    "regex_fwd": r"(TTTAGGG){2,}",
    "regex_rev": r"(CCCTAAA){2,}",
}
MOTIF_FAMILIES["insect"] = {
    "forward": ["TTAGG"],
    "regex_fwd": r"(TTAGG){2,}",
    "regex_rev": r"(CCTAA){2,}",
}


def get_taxon_families(taxon="other"):
    """Return the list of motif family names to use for a given taxon preset.

    Args:
        taxon: One of 'vertebrate', 'animal', 'plant', 'insect', 'fungal', 'other'

    Returns:
        list of str: Motif family keys from MOTIF_FAMILIES
    """
    preset = TAXON_PRESETS.get(taxon, TAXON_PRESETS["other"])
    return [f for f in preset["primary_families"] if f in MOTIF_FAMILIES]


# ── Core utility functions ───────────────────────────────────────────────────

def revcomp(seq):
    """Return reverse complement of a DNA sequence."""
    table = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(table)[::-1]


def read_fasta(path):
    """Parse FASTA and return dict {name: uppercase_sequence}."""
    seqs = {}
    name = None
    buf = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip().upper())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def write_fasta(seqs_dict, path, wrap=60):
    """Write dict {name: seq} to FASTA file."""
    with open(path, "w") as f:
        for name, seq in seqs_dict.items():
            f.write(f">{name}\n")
            for i in range(0, len(seq), wrap):
                f.write(seq[i:i + wrap] + "\n")


def extract_ends(seq, window):
    """Extract terminal sequences from a contig.

    Uses min(window, len//2) to avoid overlap of left and right ends.
    Returns ("", "") if the usable window is < 10 bp.
    """
    w = min(window, len(seq) // 2)
    if w < 10:
        return ("", "")
    return (seq[:w], seq[-w:])


# ── K-mer discovery (enrichment-based, matching TACO.sh) ────────────────────

def kmer_frequencies(seq, k):
    """Count k-mer occurrences in a sequence."""
    counts = Counter()
    for i in range(len(seq) - k + 1):
        counts[seq[i:i + k]] += 1
    return counts


def canonicalize_kmer(kmer):
    """Return the canonical form of a k-mer (min of all rotations of fwd+rev)."""
    rc = revcomp(kmer)
    candidates = []
    for x in (kmer, rc):
        doubled = x + x
        for i in range(len(x)):
            candidates.append(doubled[i:i + len(x)])
    return min(candidates)


def discover_motifs_python(end_seqs, kmin=4, kmax=30, top_n=5):
    """Discover enriched telomeric motifs from contig-end sequences.

    Uses enrichment scoring (observed / expected) with canonical k-mer rotation.
    Matches TACO.sh discover_motifs_python logic exactly.

    Args:
        end_seqs: List of end-region sequence strings
        kmin: Minimum k-mer size
        kmax: Maximum k-mer size (capped at 30)
        top_n: Number of top motifs to return

    Returns:
        list of (motif_string, count, frequency) tuples
    """
    total_bases = sum(len(s) for s in end_seqs)
    if not total_bases:
        return []

    best = []
    for k in range(kmin, min(kmax + 1, 31)):
        # Count all k-mers across all end sequences
        pool = Counter()
        for s in end_seqs:
            pool.update(kmer_frequencies(s, k))

        # Canonicalize k-mers
        canonical_counts = Counter()
        for km, c in pool.items():
            canonical_counts[canonicalize_kmer(km)] += c

        # Expected count under uniform distribution
        expected = total_bases / (4 ** k) if 4 ** k > 0 else 1

        # Score by enrichment
        for ck, c in canonical_counts.most_common(20):
            if c < 5:
                continue
            enrichment = c / max(expected, 0.01)
            if enrichment > 5:
                freq = c / max(total_bases - k + 1, 1)
                best.append((ck, c, freq, enrichment, k))

    # Sort by enrichment descending, deduplicate
    best.sort(key=lambda x: x[3], reverse=True)
    seen = set()
    result = []
    for ck, c, freq, enrichment, k in best:
        if ck in seen:
            continue
        seen.add(ck)
        result.append((ck, c, freq))
        if len(result) >= top_n:
            break
    return result


def build_regex_for_motif(motif):
    """Build forward and reverse complement regex requiring >= 2 tandem repeats.

    Args:
        motif: DNA motif string (literal, not a regex pattern)

    Returns:
        tuple: (fwd_regex_string, rev_regex_string)
    """
    rc = revcomp(motif)
    fwd = "(" + re.escape(motif) + "){2,}"
    rev = "(" + re.escape(rc) + "){2,}"
    return fwd, rev


# ── End-region scoring ───────────────────────────────────────────────────────

def score_end(end_seq, score_window, user_motif_patterns, family_patterns):
    """Score a terminal sequence for telomeric content.

    Matches TACO.sh score_end() exactly:
      - Tests all user motif patterns and all MOTIF_FAMILIES patterns
      - Selects the pattern with highest covered_bp
      - Computes composite score with proper normalization

    Args:
        end_seq: Terminal sequence to score (should be score_window-sized)
        score_window: Window size for scoring
        user_motif_patterns: list of (label, (fwd_regex_str, rev_regex_str)) tuples
        family_patterns: dict matching MOTIF_FAMILIES structure

    Returns:
        dict with keys: covered_bp, longest_run, density, best_family,
                        distance_to_end, raw_score
    """
    if not end_seq:
        return {
            "covered_bp": 0, "longest_run": 0, "density": 0.0,
            "best_family": "none", "distance_to_end": -1, "raw_score": 0.0,
        }

    w = min(score_window, len(end_seq))
    ws = end_seq[:w]

    best = {
        "covered_bp": 0, "longest_run": 0, "density": 0.0,
        "best_family": "none", "distance_to_end": w, "raw_score": 0.0,
    }

    # Collect all patterns to test
    all_patterns = []
    for label, (fwd, rev) in user_motif_patterns:
        all_patterns.append((label, fwd))
        all_patterns.append((label, rev))
    for fam_name, fam_data in family_patterns.items():
        all_patterns.append((fam_name, fam_data["regex_fwd"]))
        all_patterns.append((fam_name, fam_data["regex_rev"]))

    for label, pattern in all_patterns:
        try:
            covered = 0
            longest = 0
            min_start = w
            for m in re.finditer(pattern, ws, re.IGNORECASE):
                span = m.end() - m.start()
                covered += span
                longest = max(longest, span)
                min_start = min(min_start, m.start())
            density = covered / max(w, 1)
            if covered > best["covered_bp"]:
                best.update(
                    covered_bp=covered,
                    longest_run=longest,
                    density=density,
                    best_family=label,
                    distance_to_end=min_start,
                )
        except re.error:
            continue

    # Composite score with TACO.sh normalizations
    nl = min(best["longest_run"] / max(w * 0.3, 1), 1.0)
    nd = max(1.0 - best["distance_to_end"] / max(w, 1), 0.0)
    nc = min(best["covered_bp"] / max(w * 0.5, 1), 1.0)

    best["raw_score"] = (
        0.40 * best["density"] +
        0.30 * nl +
        0.20 * nd +
        0.10 * nc
    )
    return best


def classify_contig(left_score, right_score, strong_thr=0.25, weak_thr=0.08):
    """Classify contig based on telomere scores.

    Args:
        left_score: Raw score dict from score_end (left terminus)
        right_score: Raw score dict from score_end (right terminus)
        strong_thr: Threshold for strict_t2t / single_tel_strong
        weak_thr: Threshold for telomere_supported

    Returns:
        str: Classification tier
    """
    ls = left_score["raw_score"]
    rs = right_score["raw_score"]
    if ls >= strong_thr and rs >= strong_thr:
        return "strict_t2t"
    elif ls >= strong_thr or rs >= strong_thr:
        return "single_tel_strong"
    elif ls >= weak_thr or rs >= weak_thr:
        return "telomere_supported"
    return "none"


# ── Main detection function ──────────────────────────────────────────────────

def detect_telomeres(fasta_path, mode="hybrid", user_motif=None, end_window=5000,
                     score_window=500, kmer_min=4, kmer_max=30, threads=1,
                     taxon="other"):
    """Main hybrid telomere detection matching TACO.sh logic.

    Steps:
      1. Read FASTA, uppercase all sequences
      2. Extract contig ends at end_window size for k-mer discovery
      3. Extract contig ends at score_window size for scoring
      4. Discover motifs via tidk (if available) or enrichment-based python discovery
      5. Score each contig's left and right terminus independently
      6. Classify into strict_t2t / single_tel_strong / telomere_supported / none

    Args:
        fasta_path: Path to input FASTA
        mode: "known", "auto", or "hybrid"
        user_motif: Optional user-provided motif string
        end_window: Terminal window for k-mer discovery (default 5000)
        score_window: Scoring window size (default 500)
        kmer_min: Minimum k-mer size for discovery
        kmer_max: Maximum k-mer size for discovery
        threads: Number of threads for tidk
        taxon: Taxon preset for motif family selection (default "other")

    Returns:
        list of dicts: {contig, length, left_score, right_score, classification}
    """
    if not os.path.isfile(fasta_path) or os.path.getsize(fasta_path) == 0:
        return []

    contigs = read_fasta(fasta_path)
    if not contigs:
        return []

    # Extract ends for discovery (end_window) and scoring (score_window)
    end_left_disc = []
    end_right_disc = []
    contig_score_ends = {}  # name -> (left_score_end, right_score_end)
    contig_lengths = {}

    for name, seq in contigs.items():
        contig_lengths[name] = len(seq)
        ld, rd = extract_ends(seq, end_window)
        end_left_disc.append(ld)
        end_right_disc.append(rd)
        ls, rs = extract_ends(seq, score_window)
        contig_score_ends[name] = (ls, rs)

    all_ends = end_left_disc + end_right_disc

    # ── Build motif pattern list ──
    user_motif_patterns = []  # list of (label, (fwd_regex, rev_regex))
    discovered_motifs = []

    if user_motif and mode in ("known", "hybrid"):
        fwd, rev = build_regex_for_motif(user_motif)
        user_motif_patterns.append(("user:" + user_motif, (fwd, rev)))

    if mode in ("auto", "hybrid"):
        # Try tidk first
        tidk_path = shutil.which("tidk")
        if tidk_path:
            tf = tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False)
            try:
                for i, s in enumerate(all_ends):
                    if s:
                        tf.write(f">end_{i}\n{s}\n")
                tf.close()
                try:
                    proc = subprocess.run(
                        [tidk_path, "explore",
                         "--minimum", str(kmer_min),
                         "--maximum", str(min(kmer_max, 15)),
                         "--fasta", tf.name],
                        capture_output=True, text=True, timeout=300,
                    )
                    if proc.returncode == 0 and proc.stdout.strip():
                        for line in proc.stdout.strip().split("\n"):
                            parts = line.strip().split("\t")
                            if len(parts) >= 2 and not parts[0].startswith("#"):
                                try:
                                    discovered_motifs.append(
                                        (parts[0].strip(), int(parts[1]), 0.0))
                                except (ValueError, IndexError):
                                    continue
                except Exception:
                    pass
            finally:
                os.unlink(tf.name)

        # Fallback to Python enrichment-based discovery
        if not discovered_motifs:
            discovered_motifs = discover_motifs_python(
                all_ends, kmin=kmer_min, kmax=kmer_max, top_n=5)

        # Build regex patterns for discovered motifs
        for motif_str, cnt, freq in discovered_motifs:
            fwd, rev = build_regex_for_motif(motif_str)
            user_motif_patterns.append(("auto:" + motif_str, (fwd, rev)))

    # Use taxon-aware subset of MOTIF_FAMILIES as family_patterns
    # When user provides --motif + --telomere-mode known, family patterns
    # are still loaded but scoring is dominated by the user pattern.
    # For taxon-aware mode, only load families relevant to the taxon.
    taxon_fams = get_taxon_families(taxon)
    if taxon_fams:
        family_patterns = {k: v for k, v in MOTIF_FAMILIES.items() if k in taxon_fams}
    else:
        family_patterns = dict(MOTIF_FAMILIES)

    # ── Score and classify each contig ──
    results = []
    for name, seq in contigs.items():
        left_end, right_end = contig_score_ends[name]
        left_result = score_end(left_end, score_window,
                                user_motif_patterns, family_patterns)
        right_result = score_end(right_end, score_window,
                                 user_motif_patterns, family_patterns)
        tier = classify_contig(left_result, right_result)
        results.append({
            "contig": name,
            "length": contig_lengths[name],
            "left_score": round(left_result["raw_score"], 4),
            "right_score": round(right_result["raw_score"], 4),
            "classification": tier,
        })

    return results


# ── Output writers ───────────────────────────────────────────────────────────

def write_detection_outputs(results, fasta_path, out_prefix):
    """Write telomere detection output files.

    Writes:
      - {out_prefix}.telomere_end_scores.tsv — per-contig scores and tier
      - {out_prefix}.telo_metrics.tsv — tier count summary
      - {out_prefix}.telo.fasta — contigs with any telomeric signal
      - {out_prefix}.telo.list — backward-compatible telo list

    Args:
        results: List of detection result dicts from detect_telomeres()
        fasta_path: Path to input FASTA (for extracting sequences)
        out_prefix: Output file prefix
    """
    seqs = read_fasta(fasta_path)

    # Write .telomere_end_scores.tsv (column name "tier" for compatibility)
    scores_path = f"{out_prefix}.telomere_end_scores.tsv"
    with open(scores_path, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["contig", "length", "left_score", "right_score", "tier"],
            delimiter="\t",
        )
        w.writeheader()
        for r in results:
            w.writerow({
                "contig": r["contig"],
                "length": r.get("length", len(seqs.get(r["contig"], ""))),
                "left_score": r["left_score"],
                "right_score": r["right_score"],
                "tier": r["classification"],
            })

    # Write .telo_metrics.tsv
    counts = Counter(r["classification"] for r in results)
    metrics_path = f"{out_prefix}.telo_metrics.tsv"
    with open(metrics_path, "w") as f:
        f.write(f"strict_t2t\t{counts['strict_t2t']}\n")
        f.write(f"single_tel_strong\t{counts['single_tel_strong']}\n")
        f.write(f"telomere_supported\t{counts['telomere_supported']}\n")
        f.write(f"total_telomeric\t{sum(1 for r in results if r['classification'] != 'none')}\n")

    # Write .telo.fasta (contigs with any telomeric signal)
    telo_names = set(r["contig"] for r in results if r["classification"] != "none")
    telo_fasta_path = f"{out_prefix}.telo.fasta"
    with open(telo_fasta_path, "w") as f:
        for name, seq_str in seqs.items():
            if name in telo_names:
                f.write(f">{name}\n")
                for i in range(0, len(seq_str), 60):
                    f.write(seq_str[i:i + 60] + "\n")

    # Write backward-compatible .telo.list
    weak_thr = THRESHOLDS["telomere_supported"]
    telo_list_path = f"{out_prefix}.telo.list"
    with open(telo_list_path, "w") as f:
        for r in results:
            if r["classification"] != "none":
                length = r.get("length", len(seqs.get(r["contig"], "")))
                lp = 0 if r["left_score"] >= weak_thr else length
                rp = length if r["right_score"] >= weak_thr else 0
                f.write(f"{r['contig']}\t{lp}\t{rp}\t{length}\n")

    print(f"[INFO] Classification: strict_t2t={counts['strict_t2t']} "
          f"single_tel_strong={counts['single_tel_strong']} "
          f"telomere_supported={counts['telomere_supported']}")


# ── CLI entry point ──────────────────────────────────────────────────────────

def main():
    """Command-line interface for telomere detection."""
    parser = argparse.ArgumentParser(description="Hybrid telomere detection")
    parser.add_argument("--fasta", required=True, help="Input FASTA file")
    parser.add_argument("--mode", default="hybrid", choices=["known", "auto", "hybrid"])
    parser.add_argument("--out-prefix", required=True, help="Output file prefix")
    parser.add_argument("--motif", default="", help="User-provided motif")
    parser.add_argument("--end-window", type=int, default=5000)
    parser.add_argument("--score-window", type=int, default=500)
    parser.add_argument("--kmer-min", type=int, default=4)
    parser.add_argument("--kmer-max", type=int, default=30)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--strong-threshold", type=float, default=0.25)
    parser.add_argument("--weak-threshold", type=float, default=0.08)

    args = parser.parse_args()

    if not os.path.isfile(args.fasta) or os.path.getsize(args.fasta) == 0:
        for suffix in [".telomere_end_scores.tsv", ".telo.list",
                       ".telo.fasta", ".telo_metrics.tsv"]:
            Path(args.out_prefix + suffix).touch()
        sys.exit(0)

    results = detect_telomeres(
        args.fasta,
        mode=args.mode,
        user_motif=args.motif if args.motif else None,
        end_window=args.end_window,
        score_window=args.score_window,
        kmer_min=args.kmer_min,
        kmer_max=args.kmer_max,
        threads=args.threads,
    )

    if not results:
        for suffix in [".telomere_end_scores.tsv", ".telo.list",
                       ".telo.fasta", ".telo_metrics.tsv"]:
            Path(args.out_prefix + suffix).touch()
        sys.exit(0)

    write_detection_outputs(results, args.fasta, args.out_prefix)


if __name__ == "__main__":
    main()
