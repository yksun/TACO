"""Backbone selection and scoring for TACO."""
import os
import sys
import argparse
import csv

from taco.utils import ALL_ASSEMBLERS
from taco.policy import DEFAULT_MODE, VALID_MODES, score_assembly
ASSEMBLERS = ALL_ASSEMBLERS

#: Retained for backward compatibility with callers that introspect it.  The
#: authoritative weights now live in taco.policy.TAXON_OVERRIDES.
TAXON_WEIGHTS = {
    "fungal": {"dup_penalty": 600, "t2t_weight": 350, "n50_weight": 150, "frag_penalty": 30},
    "plant":      {"dup_penalty": 300, "t2t_weight": 200, "n50_weight": 150, "frag_penalty": 50},
    "vertebrate": {"dup_penalty": 500, "t2t_weight": 200, "n50_weight": 200, "frag_penalty": 40},
    "animal":     {"dup_penalty": 500, "t2t_weight": 200, "n50_weight": 200, "frag_penalty": 40},
    "insect":     {"dup_penalty": 500, "t2t_weight": 300, "n50_weight": 150, "frag_penalty": 30},
    "other":      {"dup_penalty": 500, "t2t_weight": 300, "n50_weight": 150, "frag_penalty": 30},
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


def _compute_score(metrics, taxon="other", genomesize=None,
                   assembly_mode=DEFAULT_MODE):
    """Compute assembly score under the active representation profile.

    Delegates to :func:`taco.policy.score_assembly` so this path and the one in
    :mod:`taco.steps` cannot disagree about what a good assembly is.  The
    metric names here are this module's CSV names; they are translated to the
    policy's names below.

    Args:
        metrics: dict of {metric: value} for a single assembly
        taxon: Taxonomy (fungal, plant, vertebrate, animal, insect, other)
        genomesize: Expected haploid genome size in bp (optional); used for the
            size term
        assembly_mode: 'primary' (nonredundant) or 'full' (alternate sequence
            retained)

    Returns:
        float: Computed score
    """
    # `contigs` and `N50` default to 1 rather than 0 to preserve the historical
    # behaviour of this entry point on rows that omit them.
    score, _ = score_assembly(
        {
            "busco_s": metrics.get("BUSCO_S", 0),
            "busco_d": metrics.get("BUSCO_D", 0),
            "busco_c": metrics.get("BUSCO_C", 0),
            "t2t": metrics.get("T2T", 0),
            "single_tel": metrics.get("single", 0),
            "merqury_comp": metrics.get("MerquryComp", 0),
            "merqury_qv": metrics.get("MerquryQV", 0),
            "contigs": metrics.get("contigs", 1),
            "n50": metrics.get("N50", 1),
            "total_len": metrics.get("total_length", 0),
        },
        mode=assembly_mode, taxon=taxon,
        expected_haploid=genomesize or 0,
    )
    return score


def select_backbone(info, mode="smart", taxon="other", genomesize=None,
                    assembly_mode=DEFAULT_MODE):
    """Select best backbone assembler using smart scoring with taxon-aware weights.

    Args:
        info: dict of {assembler: {metric: value}}
        mode: Selection mode ('smart' currently)
        taxon: Taxonomy preset (fungal, plant, vertebrate, animal, insect, other)
        genomesize: Expected haploid genome size in bp (optional)
        assembly_mode: 'primary' (default, nonredundant representation) or
            'full' (alternate/haplotype sequence retained)

    Returns:
        str: Selected assembler name
    """
    if not info:
        return None

    best_assembler = None
    best_score = float('-inf')

    for assembler, metrics in info.items():
        score = _compute_score(metrics, taxon=taxon, genomesize=genomesize,
                               assembly_mode=assembly_mode)

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
    parser.add_argument("--assembly-mode", dest="assembly_mode",
                        default=DEFAULT_MODE, choices=list(VALID_MODES),
                        help="Assembly representation to select for: 'primary' "
                             "(default, nonredundant) or 'full' (alternate/"
                             "haplotype sequence retained)")

    args = parser.parse_args()

    # Parse assembly info
    info = parse_assembly_info(args.info_csv)

    # Select backbone
    selected = select_backbone(info, mode=args.mode,
                               assembly_mode=args.assembly_mode)

    # Write debug output if requested
    if args.debug_tsv:
        with open(args.debug_tsv, 'w') as f:
            f.write("assembler\tscore\n")
            for assembler, metrics in info.items():
                score = _compute_score(metrics, assembly_mode=args.assembly_mode)
                f.write(f"{assembler}\t{score:.2f}\n")

    # Write decision output if requested
    if args.decision_txt:
        with open(args.decision_txt, 'w') as f:
            f.write(f"Selected backbone assembler: {selected}\n")
            f.write(f"Assembly representation mode: {args.assembly_mode}\n")

    # Print to stdout
    print(selected)


if __name__ == "__main__":
    main()
