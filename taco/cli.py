"""Command-line argument parsing for TACO."""
import argparse
import sys


STEP_NAMES = {
    0: "Input QC and validation",
    1: "HiCanu assembly",
    2: "NextDenovo assembly",
    3: "Peregrine assembly",
    4: "IPA assembly",
    5: "Flye assembly",
    6: "Hifiasm assembly",
    7: "LJA assembly",
    8: "MBG assembly",
    9: "Raven assembly",
    10: "Normalize + QC comparison (BUSCO + Telomere + QUAST + Merqury)",
    11: "Build telomere pool (quickmerge + validation)",
    12: "Backbone selection and refinement",
    13: "Final QC (BUSCO + Telomere + QUAST + Merqury on final)",
    14: "Final comparison report and cleanup",
    15: "Assembly-only comparison and cleanup",
}


def expand_steps(step_str):
    """Expand step string like '1,3-5,9-11' into sorted list of integers."""
    if not step_str:
        return []
    steps = set()
    for part in step_str.split(','):
        part = part.strip()
        if '-' in part:
            try:
                s, e = map(int, part.split('-'))
                if s > e: raise ValueError(f"Invalid range: {part}")
                steps.update(range(s, e + 1))
            except ValueError as ex:
                raise ValueError(f"Invalid step range: {part}") from ex
        else:
            try:
                steps.add(int(part))
            except ValueError:
                raise ValueError(f"Invalid step number: {part}")
    return sorted(list(steps))


TAXON_BUSCO_LINEAGE = {
    "fungal": "ascomycota_odb10",
    "plant": "embryophyta_odb10",
    "vertebrate": "vertebrata_odb10",
    "animal": "metazoa_odb10",
    "insect": "insecta_odb10",
    "other": None,  # require explicit --busco or warn
}


def parse_args():
    """Parse command-line arguments for TACO."""
    parser = argparse.ArgumentParser(prog='TACO', description='TACO v1.3.1 - Telomere-Aware Contig Optimization')
    parser.add_argument('-g', '--genomesize', type=str, required=True, help='Estimated genome size')
    parser.add_argument('-t', '--threads', type=int, required=True, help='Number of threads')
    parser.add_argument('--fastq', type=str, required=True, help='Path to input FASTQ')
    parser.add_argument('-m', '--motif', type=str,
                        help='Telomere motif override (optional; usually not needed — use --taxon instead)')
    parser.add_argument('--taxon', choices=['vertebrate', 'animal', 'plant', 'insect', 'fungal', 'other'],
                        default='other', help='Taxonomy preset for telomere motif priors (default: other)')
    parser.add_argument('--platform', choices=['pacbio-hifi', 'nanopore', 'pacbio'], default='pacbio-hifi')
    parser.add_argument('--reference', '-ref', type=str, help='Reference FASTA for comparison')
    parser.add_argument('-s', '--steps', type=str, help='Steps to run')
    parser.add_argument('--assembly-only', action='store_true')
    parser.add_argument('--telomere-mode', choices=['known', 'auto', 'hybrid'], default='hybrid')
    parser.add_argument('--telo-end-window', type=int, default=5000)
    parser.add_argument('--telo-score-window', type=int, default=500)
    parser.add_argument('--telo-kmer-min', type=int, default=4)
    parser.add_argument('--telo-kmer-max', type=int, default=30)
    parser.add_argument('--auto-mode', choices=['smart', 'n50'], default='smart')
    parser.add_argument('--choose', nargs='?', const='__prompt__')
    parser.add_argument('--busco', type=str, help='BUSCO lineage database (e.g., ascomycota_odb10); if not provided, defaults based on --taxon')
    parser.add_argument('--merqury', action='store_true',
                        help='Force-enable Merqury. If no --merqury-db is provided, TACO builds a reads .meryl database when meryl is installed.')
    parser.add_argument('--merqury-db', type=str,
                        help='Path to an existing reads .meryl database for Merqury')
    parser.add_argument('--merqury-k', type=str, default="auto",
                        help='K-mer size for Merqury database (default: auto; '
                             'uses Merqury best_k.sh or a genome-size fallback)')
    parser.add_argument('--no-merqury', action='store_true')
    parser.add_argument('--no-purge-dups', action='store_true', help='Skip purge_dups after refinement')
    parser.add_argument('--no-polish', action='store_true', help='Skip automatic polishing after refinement')
    parser.add_argument('--no-coverage-qc', action='store_true', help='Skip final coverage QC')
    parser.add_argument('--benchmark', action='store_true',
                        help='Write optional step timing/provenance files to benchmark_logs/')
    parser.add_argument('--allow-t2t-replace', action='store_true',
                        help='Allow rescue donors to replace immutable Tier 1 (protected T2T) contigs. '
                             'Disabled by default for safety. Use only if you have strong reason to '
                             'believe a donor is a better T2T contig than the existing one.')
    parser.add_argument('--version', action='version', version='TACO v1.3.1')
    
    args = parser.parse_args()

    # Set default BUSCO lineage based on taxon if not explicitly provided
    if args.busco is None:
        default_lineage = TAXON_BUSCO_LINEAGE.get(args.taxon)
        if default_lineage:
            args.busco = default_lineage
        # if "other" taxon with no --busco, leave as None (pipeline will warn)

    if args.telomere_mode == 'known' and not args.motif:
        parser.error("--telomere-mode known requires --motif")
    if args.no_merqury:
        args.merqury = False
    elif args.merqury or args.merqury_db:
        args.merqury = True
    
    if args.assembly_only:
        # Steps 0-10, 15: assemblers + normalize/QC + assembly-only report
        args.steps = list(range(0, 11)) + [15]
    elif args.steps:
        try:
            args.steps = expand_steps(args.steps)
        except ValueError as e:
            parser.error(str(e))
    else:
        # Full mode: Steps 0-14
        args.steps = list(range(0, 15))
    
    for s in args.steps:
        if s < 0 or s > 15:
            parser.error(
                f"Invalid step: {s}. TACO v1.3.0 uses steps 0-15. "
                f"Full mode: 0-14. Assembly-only: 0-10, 15. "
                f"Resume refinement: -s 12-14.")
    return args
