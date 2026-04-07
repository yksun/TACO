"""Backbone selection and scoring for TACO."""
import os
import sys
import argparse
import csv
import math

ASSEMBLERS = ["canu", "external", "flye", "ipa", "nextDenovo", "peregrine", "hifiasm"]


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


def select_backbone(info, mode="smart"):
    """Select best backbone assembler using smart scoring.

    Args:
        info: dict of {assembler: {metric: value}}
        mode: Selection mode ('smart' currently)

    Returns:
        str: Selected assembler name
    """
    if not info:
        return None

    best_assembler = None
    best_score = float('-inf')

    for assembler, metrics in info.items():
        # Extract metrics with defaults
        busco_s = float(metrics.get("BUSCO_S", 0))
        t2t = float(metrics.get("T2T", 0))
        single = float(metrics.get("single", 0))
        merqury_comp = float(metrics.get("MerquryComp", 0))
        merqury_qv = float(metrics.get("MerquryQV", 0))
        contigs = float(metrics.get("contigs", 1))
        n50 = float(metrics.get("N50", 1))

        # Smart scoring formula
        score = (
            busco_s * 1000 +
            t2t * 300 +
            single * 150 +
            merqury_comp * 200 +
            merqury_qv * 20 -
            contigs * 30 +
            (math.log10(n50) if n50 > 0 else 0) * 150
        )

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
                busco_s = float(metrics.get("BUSCO_S", 0))
                t2t = float(metrics.get("T2T", 0))
                single = float(metrics.get("single", 0))
                merqury_comp = float(metrics.get("MerquryComp", 0))
                merqury_qv = float(metrics.get("MerquryQV", 0))
                contigs = float(metrics.get("contigs", 1))
                n50 = float(metrics.get("N50", 1))

                score = (
                    busco_s * 1000 +
                    t2t * 300 +
                    single * 150 +
                    merqury_comp * 200 +
                    merqury_qv * 20 -
                    contigs * 30 +
                    (math.log10(n50) if n50 > 0 else 0) * 150
                )
                f.write(f"{assembler}\t{score:.2f}\n")

    # Write decision output if requested
    if args.decision_txt:
        with open(args.decision_txt, 'w') as f:
            f.write(f"Selected backbone assembler: {selected}\n")

    # Print to stdout
    print(selected)


if __name__ == "__main__":
    main()
