"""Command-line argument parsing for TACO."""
import argparse
import sys


STEP_NAMES = {
    1: "HiCanu assembly",
    2: "NextDenovo assembly",
    3: "Peregrine assembly",
    4: "IPA assembly",
    5: "Flye assembly",
    6: "Hifiasm assembly",
    7: "Copy and normalize all assemblies",
    8: "BUSCO on all assemblies",
    9: "Telomere detection and scoring",
    10: "Build optimized telomere pool",
    11: "QUAST for assembler comparison",
    12: "Backbone selection and refinement",
    13: "BUSCO on final assembly",
    14: "Telomere analysis of final assembly",
    15: "QUAST on final assembly",
    16: "Final comparison report",
    17: "Cleanup temporary files",
    18: "Assembly-only comparison summary",
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


def parse_args():
    """Parse command-line arguments for TACO."""
    parser = argparse.ArgumentParser(prog='TACO', description='TACO v1.0.0 - Telomere-Aware Contig Optimization')
    parser.add_argument('-g', '--genomesize', type=str, required=True, help='Estimated genome size')
    parser.add_argument('-t', '--threads', type=int, required=True, help='Number of threads')
    parser.add_argument('--fastq', type=str, required=True, help='Path to input FASTQ')
    parser.add_argument('-m', '--motif', type=str, help='Telomere motif')
    parser.add_argument('--platform', choices=['pacbio-hifi', 'nanopore', 'pacbio'], default='pacbio-hifi')
    parser.add_argument('--fasta', type=str, help='External FASTA')
    parser.add_argument('-s', '--steps', type=str, help='Steps to run')
    parser.add_argument('--assembly-only', action='store_true')
    parser.add_argument('--telomere-mode', choices=['known', 'auto', 'hybrid'], default='hybrid')
    parser.add_argument('--telo-end-window', type=int, default=5000)
    parser.add_argument('--telo-score-window', type=int, default=500)
    parser.add_argument('--telo-kmer-min', type=int, default=4)
    parser.add_argument('--telo-kmer-max', type=int, default=30)
    parser.add_argument('--auto-mode', choices=['smart', 'n50'], default='smart')
    parser.add_argument('--choose', nargs='?', const='__prompt__')
    parser.add_argument('--busco', nargs='?', const='ascomycota_odb10')
    parser.add_argument('--merqury', action='store_true')
    parser.add_argument('--merqury-db', type=str)
    parser.add_argument('--no-merqury', action='store_true')
    parser.add_argument('--version', action='version', version='TACO v1.0.0')
    
    args = parser.parse_args()
    
    if args.telomere_mode == 'known' and not args.motif:
        parser.error("--telomere-mode known requires --motif")
    if args.no_merqury:
        args.merqury = False
    elif args.merqury or args.merqury_db:
        args.merqury = True
    
    if args.assembly_only:
        args.steps = [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 18]
    elif args.steps:
        try:
            args.steps = expand_steps(args.steps)
        except ValueError as e:
            parser.error(str(e))
    else:
        args.steps = list(range(1, 18))
    
    for s in args.steps:
        if s < 1 or s > 18:
            parser.error(f"Invalid step: {s}")
    return args
