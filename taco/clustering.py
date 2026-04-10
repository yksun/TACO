"""Minimap2-based contig clustering with UnionFind for TACO."""
import os
import sys
import argparse
import csv

try:
    from taco.utils import read_fasta, write_fasta
except ImportError:
    try:
        from utils import read_fasta, write_fasta
    except ImportError:
        pass


class UnionFind:
    """Union-Find data structure for clustering."""

    def __init__(self):
        """Initialize empty Union-Find."""
        self.parent = {}
        self.rank = {}

    def find(self, x):
        """Find the root parent of x with path compression."""
        if x not in self.parent:
            self.parent[x] = x
            self.rank[x] = 0
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, x, y):
        """Union two elements by rank."""
        rx, ry = self.find(x), self.find(y)
        if rx == ry:
            return
        if self.rank[rx] < self.rank[ry]:
            rx, ry = ry, rx
        self.parent[ry] = rx
        if self.rank[rx] == self.rank[ry]:
            self.rank[rx] += 1


def parse_paf_and_cluster(paf_path, seq_names, min_identity=0.95, min_coverage=0.85):
    """Parse PAF file and cluster contigs based on alignment quality.

    Uses QUERY coverage only to match TACO.sh logic — a shorter contig
    fully covered by a longer one triggers union, but not the reverse.

    Args:
        paf_path: Path to PAF file
        seq_names: Set/list of valid sequence names (ensures singletons are included)
        min_identity: Minimum identity threshold (default 0.95)
        min_coverage: Minimum query coverage threshold (default 0.85)

    Returns:
        dict: {cluster_root: [member_names]}
    """
    uf = UnionFind()

    # Initialise every known sequence so singletons appear as their own cluster
    for n in seq_names:
        uf.find(n)

    # Parse PAF file
    if os.path.isfile(paf_path):
        with open(paf_path, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 11:
                    continue

                query_name = parts[0]
                query_len = int(parts[1])
                query_start = int(parts[2])
                query_end = int(parts[3])

                target_name = parts[5]

                num_matches = int(parts[9])
                alignment_len = int(parts[10])

                # Skip self-matches and zero-length alignments
                if query_name == target_name or alignment_len == 0 or query_len == 0:
                    continue
                # Skip names not in the input set
                if query_name not in uf.parent or target_name not in uf.parent:
                    continue

                identity = num_matches / alignment_len
                query_coverage = (query_end - query_start) / query_len

                if identity >= min_identity and query_coverage >= min_coverage:
                    uf.union(query_name, target_name)

    # Build clusters
    clusters = {}
    for name in seq_names:
        root = uf.find(name)
        clusters.setdefault(root, []).append(name)

    return clusters


def cluster_and_select(fasta_path, paf_path, out_fasta, out_tsv, label="",
                       min_identity=0.95, min_coverage=0.85, scores=None):
    """Cluster contigs and select the best representative per cluster.

    Selection priority:
        1. Telomere score (if ``scores`` dict provided) — highest score wins
        2. Sequence length — longest contig wins (tiebreaker)
        3. Name (alphabetical) — deterministic final tiebreaker

    Args:
        fasta_path: Path to input FASTA
        paf_path: Path to PAF alignment file
        out_fasta: Output FASTA path
        out_tsv: Output summary TSV path
        label: Label for cluster summary
        min_identity: Minimum identity threshold
        min_coverage: Minimum query coverage threshold
        scores: Optional dict {contig_name: float} of telomere scores for
                selecting the best representative.  When provided, the contig
                with the highest score is preferred; ties are broken by length.
    """
    # Read sequences and lengths
    seqs = read_fasta(fasta_path)
    lengths = {name: len(seq) for name, seq in seqs.items()}

    # Cluster using query-coverage only (matching TACO.sh)
    clusters = parse_paf_and_cluster(paf_path, set(seqs.keys()),
                                     min_identity, min_coverage)

    # Select best per cluster: (telomere_score, length, name) — highest wins
    selected = {}
    summary = []

    for root, members in sorted(clusters.items()):
        if scores:
            # Sort by (telomere_score DESC, length DESC, name ASC for determinism)
            members_sorted = sorted(
                members,
                key=lambda x: (scores.get(x, 0.0), lengths.get(x, 0), x),
                reverse=True,
            )
        else:
            # Fallback: longest contig (TACO.sh default)
            members_sorted = sorted(
                members,
                key=lambda x: (lengths.get(x, 0), x),
                reverse=True,
            )
        rep = members_sorted[0]
        selected[rep] = seqs[rep]

        summary.append({
            "cluster_id": root,
            "representative": rep,
            "rep_length": lengths[rep],
            "rep_score": f"{scores.get(rep, 0.0):.4f}" if scores else "",
            "cluster_size": len(members),
            "members": ",".join(sorted(members)),
            "label": label,
        })

    # Write output FASTA
    write_fasta(selected, out_fasta)

    # Write summary TSV
    with open(out_tsv, "w", newline="") as f:
        fieldnames = ["cluster_id", "representative", "rep_length", "rep_score",
                      "cluster_size", "members", "label"]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)

    n_in = len(seqs)
    n_out = len(selected)
    print(f"Clustering complete: {n_in} input → {len(clusters)} clusters → {n_out} representatives")
    return selected


def main():
    """Command-line interface for clustering."""
    parser = argparse.ArgumentParser(description="Minimap2-based contig clustering")
    parser.add_argument("--fasta", required=True, help="Input FASTA file")
    parser.add_argument("--paf", required=True, help="Input PAF file from minimap2")
    parser.add_argument("--out-fasta", required=True, help="Output FASTA file")
    parser.add_argument("--out-tsv", required=True, help="Output summary TSV file")
    parser.add_argument("--label", default="", help="Label for summary")
    parser.add_argument("--min-identity", type=float, default=0.95, help="Minimum identity threshold")
    parser.add_argument("--min-coverage", type=float, default=0.85, help="Minimum coverage threshold")

    args = parser.parse_args()

    cluster_and_select(
        args.fasta,
        args.paf,
        args.out_fasta,
        args.out_tsv,
        label=args.label,
        min_identity=args.min_identity,
        min_coverage=args.min_coverage,
    )


if __name__ == "__main__":
    main()
