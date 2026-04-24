"""Backbone selection and scoring for TACO."""
import os
import sys
import argparse
import csv
import math

ASSEMBLERS = ["canu", "reference", "flye", "ipa", "nextDenovo", "peregrine", "hifiasm"]

TAXON_WEIGHTS = {
    "fungal": {"dup_penalty": 600, "t2t_weight": 350, "n50_weight": 150, "frag_penalty": 30},
    "plant": {"dup_penalty": 200, "t2t_weight": 200, "n50_weight": 100, "frag_penalty": 15},
    "vertebrate": {"dup_penalty": 400, "t2t_weight": 250, "n50_weight": 120, "frag_penalty": 20},
    "animal": {"dup_penalty": 400, "t2t_weight": 250, "n50_weight": 120, "frag_penalty": 20},
    "insect": {"dup_penalty": 500, "t2t_weight": 300, "n50_weight": 130, "frag_penalty": 25},
    "other": {"dup_penalty": 400, "t2t_weight": 300, "n50_weight": 150, "frag_penalty": 30},
}


def parse_assembly_info(csv_path):
    """Parse assemblies/assembly_info.csv into dict of {assembler: {metric: value}}.

    Args:
        csv_path: Path to assembly_info.csv file

    Returns:
        dict: {assembler_name: {metric_name: metric_value}}
    """
    info = {}

    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Assembly info file not found: {csv_path}")

    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row is None or not row:
                continue

            assembler = row.get("assembler", "").strip()
            if not assembler:
                continue

            # Create entry for this assembler
            if assembler not in info:
                info[assembler] = {}

            # Store all metrics
            for key, value in row.items():
                if key and key != "assembler":
                    try:
                        # Try to convert to float
                        info[assembler][key] = float(value)
                    except (ValueError, TypeError):
                        # Keep as string if not numeric
                        info[assembler][key] = value

    return info


def _compute_score(metrics, taxon="other", genomesize=None):
    """Compute assembly score using taxon-aware weights and BUSCO D penalty.

    Args:
        metrics: dict of {metric: value} for a single assembly
        taxon: Taxonomy (fungal, plant, vertebrate, animal, insect, other)
        genomesize: Expected genome size in bp (optional); used for size deviation penalty

    Returns:
        float: Computed score
    """
    # Get taxon-specific weights
    weights = TAXON_WEIGHTS.get(taxon, TAXON_WEIGHTS["other"])
    dup_penalty = weights["dup_penalty"]
    t2t_weight = weights["t2t_weight"]
    n50_weight = weights["n50_weight"]
    frag_penalty = weights["frag_penalty"]

    # Extract metrics with defaults
    busco_s = float(metrics.get("BUSCO_S", 0))
    busco_d = float(metrics.get("BUSCO_D", 0))
    t2t = float(metrics.get("T2T", 0))
    single = float(metrics.get("single", 0))
    merqury_comp = float(metrics.get("MerquryComp", 0))
    merqury_qv = float(metrics.get("MerquryQV", 0))
    contigs = float(metrics.get("contigs", 1))
    n50 = float(metrics.get("N50", 1))

    # Size deviation penalty
    size_penalty = 0
    if genomesize and genomesize > 0:
        total_len = float(metrics.get("total_length", 0))
        if total_len > 0:
            deviation = abs(total_len - genomesize) / genomesize
            size_penalty = deviation * 500  # 500 points penalty per 100% deviation
        else:
            size_penalty = 0
    else:
        size_penalty = 0

    # Smart scoring formula with BUSCO_D penalty and taxon-aware weights
    score = (
        busco_s * 1000
        - busco_d * dup_penalty
        + merqury_comp * 200
        + merqury_qv * 20
        + t2t * t2t_weight
        + single * 150
        + (math.log10(n50) if n50 > 0 else 0) * n50_weight
        - contigs * frag_penalty
        - size_penalty
    )

    return score


def select_backbone(info, mode="smart", taxon="other", genomesize=None):
    """Select best backbone assembler using smart scoring with taxon-aware weights.

    Args:
        info: dict of {assembler: {metric: value}}
        mode: Selection mode ('smart' currently)
        taxon: Taxonomy preset (fungal, plant, vertebrate, animal, insect, other)
        genomesize: Expected genome size in bp (optional)

    Returns:
        str: Selected assembler name
    """
    if not info:
        return None

    best_assembler = None
    best_score = float('-inf')

    for assembler, metrics in info.items():
        score = _compute_score(metrics, taxon=taxon, genomesize=genomesize)

        if score > best_score:
            best_score = score
            best_assembler = assembler

    return best_assembler


def main():
    """Command-line interface for backbone selection."""
    parser = argparse.ArgumentParser(description="Backbone assembler selection")
    parser.add_argument("--info-csv", required=True, help="Path to assembly_info.csv")
    parser.add_argument("--mode", default="smart", choices=["smart"],
                        help="Selection mode")
    parser.add_argument("--debug-tsv", help="Optional debug TSV output with scores")
    parser.add_argument("--decision-txt", help="Optional decision text output")

    args = parser.parse_args()

    # Parse assembly info
    info = parse_assembly_info(args.info_csv)

    # Select backbone
    selected = select_backbone(info, mode=args.mode)

    # Write debug output if requested
    if args.debug_tsv:
        with open(args.debug_tsv, 'w') as f:
            f.write("assembler\tscore\n")
            for assembler, metrics in info.items():
                score = _compute_score(metrics)
                f.write(f"{assembler}\t{score:.2f}\n")

    # Write decision output if requested
    if args.decision_txt:
        with open(args.decision_txt, 'w') as f:
            f.write(f"Selected backbone assembler: {selected}\n")

    # Print to stdout
    print(selected)


if __name__ == "__main__":
    main()
