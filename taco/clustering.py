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
        """Find the root parent of x with path compression.

        Args:
            x: Element to find

        Returns:
            Root parent of x
        """
        if x not in self.parent:
            self.parent[x] = x
            self.rank[x] = 0
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, x, y):
        """Union two elements by rank.

        Args:
            x: First element
            y: Second element
        """
        rx, ry = self.find(x), self.find(y)
        if rx == ry:
            return
        if self.rank[rx] < self.rank[ry]:
            rx, ry = ry, rx
        self.parent[ry] = rx
        if self.rank[rx] == self.rank[ry]:
            self.rank[rx] += 1


def read_fasta_lengths(fasta_path):
    """Read FASTA file and return dict of {name: length}.

    Args:
        fasta_path: Path to FASTA file

    Returns:
        dict: {contig_name: sequence_length}
    """
    seqs = read_fasta(fasta_path)
    return {name: len(seq) for name, seq in seqs.items()}


def parse_paf_and_cluster(paf_path, min_identity=0.95, min_coverage=0.85):
    """Parse PAF file and cluster contigs based on alignment quality.

    Args:
        paf_path: Path to PAF file
        min_identity: Minimum identity threshold (default 0.95)
        min_coverage: Minimum coverage threshold (default 0.85)

    Returns:
        dict: {cluster_root: [member_names]}
    """
    uf = UnionFind()
    lengths = {}

    # Parse PAF file
    with open(paf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 11:
                continue

            query_name = parts[0]
            query_len = int(parts[1])
            query_start = int(parts[2])
            query_end = int(parts[3])

            target_name = parts[5]
            target_len = int(parts[6])
            target_start = int(parts[7])
            target_end = int(parts[8])

            num_matches = int(parts[9])
            alignment_len = int(parts[10])

            # Store lengths
            lengths[query_name] = query_len
            lengths[target_name] = target_len

            # Calculate identity
            if alignment_len == 0:
                continue
            identity = num_matches / alignment_len

            # Calculate coverage
            query_coverage = (query_end - query_start) / query_len
            target_coverage = (target_end - target_start) / target_len
            min_cov = min(query_coverage, target_coverage)

            # Cluster if thresholds met
            if identity >= min_identity and min_cov >= min_coverage:
                uf.union(query_name, target_name)

    # Build clusters
    clusters = {}
    for name in lengths:
        root = uf.find(name)
        if root not in clusters:
            clusters[root] = []
        clusters[root].append(name)

    return clusters


def cluster_and_select(fasta_path, paf_path, out_fasta, out_tsv, label="",
                       min_identity=0.95, min_coverage=0.85):
    """Cluster contigs and select longest representative per cluster.

    Args:
        fasta_path: Path to input FASTA
        paf_path: Path to PAF alignment file
        out_fasta: Output FASTA path
        out_tsv: Output summary TSV path
        label: Label for cluster summary
        min_identity: Minimum identity threshold
        min_coverage: Minimum coverage threshold
    """
    # Read sequences and lengths
    seqs = read_fasta(fasta_path)
    lengths = {name: len(seq) for name, seq in seqs.items()}

    # Cluster
    clusters = parse_paf_and_cluster(paf_path, min_identity, min_coverage)

    # Select longest per cluster
    selected = {}
    summary = []

    for root, members in clusters.items():
        # Find longest member
        longest_name = max(members, key=lambda x: lengths.get(x, 0))
        longest_len = lengths[longest_name]

        selected[longest_name] = seqs[longest_name]

        summary.append({
            "cluster_id": root,
            "representative": longest_name,
            "rep_length": longest_len,
            "cluster_size": len(members),
            "members": ",".join(sorted(members)),
            "label": label,
        })

    # Write output FASTA
    write_fasta(selected, out_fasta)

    # Write summary TSV
    with open(out_tsv, 'w', newline='') as f:
        fieldnames = ["cluster_id", "representative", "rep_length", "cluster_size", "members", "label"]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(summary)

    print(f"Clustering complete: {len(clusters)} clusters, {len(selected)} representatives")


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
