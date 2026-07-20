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
    14: "Report + cleanup (14A: full comparison / 14B: assembly-only)",
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
    # Use the broadest reasonable lineage for each taxon by default.
    # A specific sub-lineage (e.g. ascomycota_odb10) should only be used
    # when the user explicitly passes --busco <lineage>.
    "fungal": "fungi_odb10",
    "plant": "embryophyta_odb10",
    "vertebrate": "vertebrata_odb10",
    "animal": "metazoa_odb10",
    "insect": "insecta_odb10",
    "other": None,  # require explicit --busco or warn
}


def parse_args():
    """Parse command-line arguments for TACO."""
    parser = argparse.ArgumentParser(prog='TACO', description='TACO v1.3.6 - Telomere-Aware Contig Optimization')
    parser.add_argument('-g', '--genomesize', type=str, required=True, help='Estimated genome size')
    parser.add_argument('-t', '--threads', type=int, required=True, help='Number of threads')
    parser.add_argument('--fastq', type=str, required=True, help='Path to input FASTQ')
    parser.add_argument('-m', '--motif', type=str,
                        help='Telomere motif override (optional; usually not needed — use --taxon instead)')
    parser.add_argument('--taxon', choices=['vertebrate', 'animal', 'plant', 'insect', 'fungal', 'other'],
                        default='other', help='Taxonomy preset for telomere motif priors (default: other)')
    parser.add_argument('--platform', choices=['pacbio-hifi', 'nanopore', 'pacbio'], default='pacbio-hifi')
    parser.add_argument('--reference', '-ref', type=str,
                        help='Reference FASTA. Included as the "reference" row in '
                             'comparison tables (BUSCO/Telomere/QUAST/Merqury). '
                             'Excluded from backbone selection, quickmerge candidate '
                             'pool, telomere-aware refinement, polish, and purge_dups.')
    parser.add_argument('--compare', type=str,
                        help='Compare-only FASTA (e.g. an external assembly to '
                             'compare against). Same QC as --reference, AND TACO '
                             'emits a contig-to-contig comparison report '
                             '(final_results/compare_report/) against the final '
                             'merged assembly: minimap2 asm5 alignment, per-contig '
                             '1-to-1 mapping, alignment identity, and weak_regions '
                             '(stretches of final.merged.fasta with no/low support '
                             'from the compare assembly). Never used for selection, '
                             'merging, or polishing.')
    parser.add_argument('--final-fa', dest='final_fa', type=str,
                        help='Use this FASTA as the final merged assembly for steps '
                             '13/14 instead of assemblies/final.merged.fasta. Useful '
                             'for re-running final QC against an externally produced '
                             'assembly without rebuilding from step 12.')
    parser.add_argument('-s', '--steps', type=str, help='Steps to run')
    parser.add_argument('--assembly-only', action='store_true')
    parser.add_argument('--telomere-mode', choices=['known', 'auto', 'hybrid'], default='hybrid')
    parser.add_argument('--telo-end-window', type=int, default=5000)
    parser.add_argument('--telo-score-window', type=int, default=None,
                        help='Telomere scoring window (bp). Default: taxon-aware '
                             '(fungal 300; plant/vertebrate/animal 1000; other 500). '
                             'An explicit value here is always honored.')
    parser.add_argument('--telo-kmer-min', type=int, default=4)
    parser.add_argument('--telo-kmer-max', type=int, default=30)
    parser.add_argument('--auto-mode', choices=['smart', 'n50'], default='smart')
    parser.add_argument('--choose', nargs='?', const='__prompt__')
    parser.add_argument('--busco', type=str,
                        help='BUSCO lineage database (e.g., ascomycota_odb10, fungi_odb10). '
                             'Only used as-is when explicitly provided. If omitted, TACO picks a '
                             'broad default based on --taxon (fungal -> fungi_odb10, plant -> '
                             'embryophyta_odb10, vertebrate -> vertebrata_odb10, animal -> '
                             'metazoa_odb10, insect -> insecta_odb10).')
    parser.add_argument('--busco-download-path', dest='busco_download_path', type=str, default=None,
                        help='Directory where BUSCO lineage datasets are cached '
                             '(passed as --download_path to busco). Default: BUSCO_DOWNLOAD_PATH '
                             'env var if set, otherwise ./busco_downloads relative to the run dir.')
    parser.add_argument('--busco-offline-only', dest='busco_offline_only', action='store_true',
                        help='Refuse to download BUSCO lineages over the network. If the lineage '
                             'is not already cached in --busco-download-path, BUSCO will fail '
                             'instead of falling back to an online lookup.')
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
    parser.add_argument('--version', action='version', version='TACO v1.3.6')
    
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
    
    # Steps 0-10, 14: assemblers + normalize/QC + report (14B assembly-only)
    # Skips 11 (telomere pool), 12 (refinement), 13 (final QC)
    assembly_only_steps = list(range(0, 11)) + [14]
    if args.assembly_only:
        if args.steps:
            # Honor an explicit --steps subset instead of silently discarding
            # it, but keep it within the assembly-only step set.
            try:
                requested = expand_steps(args.steps)
            except ValueError as e:
                parser.error(str(e))
            allowed = set(assembly_only_steps)
            args.steps = [s for s in requested if s in allowed]
            dropped = [s for s in requested if s not in allowed]
            if dropped:
                print(f"[warn] --assembly-only ignores step(s) {dropped}: "
                      f"assembly-only mode covers steps 0-10 and 14 only.",
                      file=sys.stderr)
            if not args.steps:
                parser.error(
                    "--assembly-only combined with --steps left no runnable "
                    "steps; assembly-only mode covers steps 0-10 and 14.")
        else:
            args.steps = assembly_only_steps
    elif args.steps:
        try:
            args.steps = expand_steps(args.steps)
        except ValueError as e:
            parser.error(str(e))
    else:
        # Full mode: Steps 0-14 (14A runs automatically)
        args.steps = list(range(0, 15))
    
    for s in args.steps:
        if s < 0 or s > 14:
            parser.error(
                f"Invalid step: {s}. TACO v1.3.6 uses steps 0-14. "
                f"Full mode: 0-14. Assembly-only: 0-10, 14. "
                f"Resume from refinement: -s 12-14.")
    return args
