"""Telomere pool classification for TACO."""
import os
import sys
import argparse

try:
    from taco.telomere_detect import detect_telomeres, MOTIF_FAMILIES
    from taco.utils import read_fasta, write_fasta
except ImportError:
    try:
        from telomere_detect import detect_telomeres, MOTIF_FAMILIES
        from utils import read_fasta, write_fasta
    except ImportError:
        detect_telomeres = None
        MOTIF_FAMILIES = {}
        read_fasta = None
        write_fasta = None


def classify_pool(fasta_path, mode="hybrid", motif=None, out_prefix="pool",
                  end_window=5000, score_window=500, kmer_min=4, kmer_max=30, threads=1):
    """Run telomere detection on pooled contigs and write classified outputs.

    Args:
        fasta_path: Path to merged FASTA file
        mode: Detection mode ('known', 'auto', 'hybrid')
        motif: User-provided motif (optional)
        out_prefix: Output file prefix
        end_window: Terminal window size
        score_window: Scoring window size
        kmer_min: Minimum k-mer size
        kmer_max: Maximum k-mer size
        threads: Number of threads
    """
    if detect_telomeres is None or read_fasta is None or write_fasta is None:
        raise ImportError("Could not import required functions from telomere_detect and utils")

    # Run telomere detection
    results = detect_telomeres(
        fasta_path,
        mode=mode,
        user_motif=motif,
        end_window=end_window,
        score_window=score_window,
        kmer_min=kmer_min,
        kmer_max=kmer_max,
        threads=threads,
    )

    # Read sequences
    seqs = read_fasta(fasta_path)

    # Classify into tiers
    t2t_contigs = {}
    single_tel_contigs = {}
    telomere_supported_contigs = {}

    for result in results:
        contig_name = result["contig"]
        classification = result["classification"]
        seq = seqs[contig_name]

        if classification == "strict_t2t":
            t2t_contigs[contig_name] = seq
            single_tel_contigs[contig_name] = seq
            telomere_supported_contigs[contig_name] = seq
        elif classification == "single_tel_strong":
            single_tel_contigs[contig_name] = seq
            telomere_supported_contigs[contig_name] = seq
        elif classification == "telomere_supported":
            telomere_supported_contigs[contig_name] = seq

    # Write t2t files
    t2t_fasta = f"{out_prefix}.t2t.fasta"
    write_fasta(t2t_contigs, t2t_fasta)
    with open(f"{out_prefix}.t2t.list", 'w') as f:
        for name in t2t_contigs:
            f.write(f"{name}\n")

    # Write single_tel files
    single_tel_fasta = f"{out_prefix}.single_tel.fasta"
    write_fasta(single_tel_contigs, single_tel_fasta)
    with open(f"{out_prefix}.single_tel.list", 'w') as f:
        for name in single_tel_contigs:
            f.write(f"{name}\n")

    # Write telomere_supported files
    telo_supported_fasta = f"{out_prefix}.telomere_supported.fasta"
    write_fasta(telomere_supported_contigs, telo_supported_fasta)
    with open(f"{out_prefix}.telomere_supported.list", 'w') as f:
        for name in telomere_supported_contigs:
            f.write(f"{name}\n")

    print(f"Pool classification complete:")
    print(f"  T2T contigs: {len(t2t_contigs)}")
    print(f"  Single telomere contigs: {len(single_tel_contigs)}")
    print(f"  Telomere-supported contigs: {len(telomere_supported_contigs)}")


def main():
    """Command-line interface for pool classification."""
    parser = argparse.ArgumentParser(description="Telomere pool classification")
    parser.add_argument("--fasta", required=True, help="Input FASTA file")
    parser.add_argument("--mode", default="hybrid", choices=["known", "auto", "hybrid"],
                        help="Detection mode")
    parser.add_argument("--out-prefix", default="pool", help="Output file prefix")
    parser.add_argument("--motif", help="User-provided motif")
    parser.add_argument("--end-window", type=int, default=5000, help="Terminal window size")
    parser.add_argument("--score-window", type=int, default=500, help="Scoring window size")
    parser.add_argument("--kmer-min", type=int, default=4, help="Minimum k-mer size")
    parser.add_argument("--kmer-max", type=int, default=30, help="Maximum k-mer size")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads")

    args = parser.parse_args()

    classify_pool(
        args.fasta,
        mode=args.mode,
        motif=args.motif,
        out_prefix=args.out_prefix,
        end_window=args.end_window,
        score_window=args.score_window,
        kmer_min=args.kmer_min,
        kmer_max=args.kmer_max,
        threads=args.threads,
    )


if __name__ == "__main__":
    main()
