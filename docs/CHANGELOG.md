# TACO Changelog

All notable changes to TACO are documented here.
Versions follow [Semantic Versioning](https://semver.org/).

---

## [1.3.3] — 2026-04-30

### Coverage-guided partial T2T replacement and upstream-aligned purge_dups

Investigation of a fungal HiFi run (Step 12, 2026-04-30) where the pool
contained three "strict T2T" contigs (`contig_118` at 207 kb covering ~20 %
of a 1.04 Mb backbone contig, plus `contig_4819` and `contig_5514` at ~2 kb
each — telomere-repeat fragments, not whole chromosomes) revealed two
behavioural gaps that 1.3.3 fixes.

- **Behavioural fix.** The strict-T2T fragments correctly cannot
  replace Mb-scale backbone contigs.  However, a substantial partial T2T
  hit (e.g. 207 kb at 99.8 % identity covering 20 % of a 1.04 Mb backbone)
  used to be silently rejected because target coverage fell below the 50 %
  full-replacement floor.  TACO now drops into a **read-coverage
  diagnostic** for these cases instead of rejecting outright.  If the
  unmatched backbone region has < 30 % of the matched-region coverage
  (configurable via `CHIMERIC_COV_RATIO`), the backbone is treated as
  chimeric/duplicated and replaced by the T2T contig; otherwise the
  backbone is kept.  Tiny (≤ taxon `PARTIAL_T2T_MIN_BP`) telomere-repeat
  fragments still cannot replace large backbone contigs.
- **Taxon-aware partial-T2T thresholds.** Coverage-guided partial
  replacement now has explicit per-taxon defaults (fungi/insect: tcov ≥
  15 %, len ≥ 50 kb, qcov ≥ 80 %, identity ≥ 90 %; plant/vertebrate/
  animal: tcov ≥ 50 %, len ≥ 200 kb, qcov ≥ 85 %, identity ≥ 92 %; other:
  tcov ≥ 30 %, len ≥ 100 kb).  All four thresholds are independently
  overridable via `PARTIAL_T2T_MIN_TCOV`, `PARTIAL_T2T_MIN_BP`,
  `PARTIAL_T2T_MIN_QCOV`, `PARTIAL_T2T_MIN_IDENT`.
- **purge_dups safety harmonised with genome-size budget.** The bp-drop
  safety check now accepts large drops when they move an over-large
  assembly closer to the expected genome size *and* stay above the
  expected-size floor (`PURGE_DUPS_MIN_EXPECTED_RATIO`).  The Step 12L
  "do no harm" comparison uses the same expected-size logic, so the
  shrink warning no longer fires when purging usefully right-sizes an
  assembly.
- **Fungal purge_dups defaults aligned with upstream.** The fungal
  preset now uses `-2 -f 0.80 -l 10000 -E 15000`, matching the
  upstream binary defaults and the recommended example pipeline in
  `dfguan/purge_dups`.  This is safer than the previous short-match
  fungal tuning (`-l 5000 -E 5000`) when the selected backbone is
  already close to the expected genome size.  `get_seqs -e` (end-only
  duplication removal) remains the default, matching the upstream
  recommendation.
- **Diagnostic logging.** Tier 2 backbone contigs without a full T2T
  upgrade now report the best partial T2T hit and flag potential
  chimeric backbones explicitly in the Step 12 log.

### Final report metric completeness (carried forward)

- Final BUSCO report rows still include `BUSCO C (count)` and
  `BUSCO M (count)` for the merged assembly.
- Merqury QV reporting still searches nested Merqury output paths and
  reports `NA` (rather than blank) when QV is missing.

---

## [1.3.2] — 2026-04-28

### Final report metric completeness

- **Fixed** final BUSCO report rows.  Step 13 now writes `BUSCO C (count)`
  and `BUSCO M (count)` for the merged assembly, matching the per-assembler
  BUSCO rows already present in `assemblies/assembly_info.csv`.
- **Improved** Merqury QV reporting.  TACO now searches nested Merqury output
  paths for `.qv` files, reports `NA` instead of a silent blank when Merqury is
  enabled but no parseable QV file exists, and logs a warning naming the
  command log and QV candidates when completeness is present but QV is missing.
- **Improved** final report clarity.  Missing final QC values are now explicit
  where they reflect a missing upstream metric rather than a CSV formatting
  gap.

---

## [1.3.1] — 2026-04-25

### Merqury QC flow hardening

- **Fixed** Merqury metric parsing.  TACO now reads QV from the fourth column
  of `.qv` files and completeness from the fifth column of
  `.completeness.stats` (preferring the `all` row), avoiding accidental use of
  k-mer count columns in assembly scoring and comparison tables.
- **Improved** Merqury database handling.  Step 10 and Step 13 now share one
  resolver that finds an existing reads `.meryl` database or builds
  `merqury/reads.k{K}.meryl` from the input FASTQ when `meryl` is installed.
- **Improved** Merqury k-mer selection.  `--merqury-k` now defaults to `auto`;
  TACO uses Merqury `best_k.sh` when available and otherwise falls back to the
  published genome-size/collision-rate formula, clamped to k=17-31 for broad
  eukaryotic assemblies.  The fallback collision rate is configurable with
  `MERQURY_COLLISION_RATE` (default 0.001).
- **Improved** final-QC logic.  Step 13 is now the only final-QC step
  (BUSCO, telomere detection, QUAST, and Merqury on the final refined
  assembly).  Step 14A builds the full report from those Step 13 metrics, and
  Step 14B is reserved for `--assembly-only`.
- **Documented** the non-high-accuracy warning for Nanopore and PacBio CLR
  reads in the README and installation guide, including cautious
  interpretation of QV/completeness from non-HiFi reads.

### Benchmark logging

- **Changed** benchmark timing/provenance output to opt-in via `--benchmark`.
  Normal pipeline runs still write per-step logs, comparison tables, and
  version metadata, but no longer create or update `benchmark_logs/` unless the
  benchmark flag is used.
- **Improved** `step_benchmark.tsv` for single-run interpretation.  The table is
  reset at run start, records `success`/`failed` status plus numeric
  `exit_code`, and `run_summary.txt` now reports current-run step counts and
  total recorded step runtime.
- **Added** publication-oriented benchmark provenance files:
  `run_manifest.json`, `software_versions.tsv`, `output_manifest.tsv`, and
  `methods_note.txt`.  These capture exact command line, code commit/dirty
  state, parameters, input file metadata, tool versions, key output files, and a
  short methods note for paper supplements.  Input SHA-256 checksums remain
  opt-in with `TACO_BENCHMARK_SHA256=1` to avoid unexpectedly hashing very large
  FASTQ files.
- **Fixed** `version.txt` handling for tools without clean version commands.
  LJA command-line parse errors are no longer recorded as versions, and
  `merge_wrapper.py` falls back to the active conda/mamba `quickmerge` package
  version when the wrapper does not expose a version flag.
- **Improved** cleanup organization.  Final report cleanup now copies stable
  results into `final_results/` while leaving source files in place for easier
  reruns/resumes, writes both `final.merged.fasta` and `final_assembly.fasta`,
  copies telomere-pool products into `telomere_pool/`, and moves large
  transient work directories into structured `temp/` subfolders.
- **Improved** resume and selected-step warnings.  Tool preflight checks are now
  based on the requested steps rather than every platform-compatible assembler,
  so resume runs such as `-s 13-14` no longer warn about missing Step 1-9
  assembler binaries.  Before each resumed step, TACO warns about missing
  upstream files and can restore common inputs from `final_results/` or
  `telomere_pool/`.
- **Fixed** `purge_dups` version parsing.  Error output such as
  `[E::hit_read] can not open PAF file version` is rejected as non-version
  text, and TACO falls back to the active conda package version when available.
- **Fixed** Raven command compatibility.  Step 9 now inspects the installed
  Raven CLI before running, uses a thread flag only when the executable
  advertises one, and skips non-assembler `raven` executables with a clearer
  `raven-assembler` installation hint.  Each attempt keeps its own stderr log
  and failed attempts print a short diagnostic tail.
- **Improved** step diagnostics.  Per-step logs now include START/END markers,
  preflight restore/missing-input messages, command exit codes, elapsed command
  times, expected assembler output paths, and short tails from tool-specific
  logs when assembler, BUSCO, QUAST, purge_dups, polishing, or Merqury commands
  fail or produce empty output.
- **Improved** MBG handling.  MBG remains included in the assembler comparison
  when installed, is now included in `taco-env.yml` via Bioconda, and TACO now
  passes MBG's required odd `-k` k-mer size (default `1501`, override with
  `TACO_MBG_K`) while treating missing MBG as an optional HiFi assembler skip
  with a clearer warning instead of a broad environment warning.
- **Changed** default step flow.  Step 10 now runs assembly normalization plus
  BUSCO, telomere, QUAST, and Merqury comparison.  Default full mode runs
  Steps 0-14.  `--assembly-only` runs Steps 0-10 and 14, where Step 14 selects
  14B only because assembly-only mode is enabled.
- **Improved** output organization.  Final cleanup and assembly-only cleanup
  now move raw assembler work directories (`hicanu/`, `flye/`, `raven_out/`,
  and peers) into `temp/assemblers/` after normalized FASTAs and comparison
  tables have been written; Step 10 can normalize from those organized
  directories during resumed runs.
- **Improved** selected-step preflight.  Standalone runs now warn about the
  specific upstream files they need: Step 10 checks for assembler outputs from
  Steps 1-9 or existing normalized FASTAs, Step 12 checks for Step 10
  comparison outputs and Step 11 telomere-pool files, Step 13 checks only for
  the final refined FASTA, and Step 14 adapts its report preflight to 14A or
  14B.  Invalid requests above Step 14 now produce guidance for the current
  0-14 step layout.
- **Fixed** Merqury result discovery.  TACO now searches exact, prefix-based,
  and nested Merqury output paths for `.qv` and `.completeness.stats` files
  before writing `assembly.merqury.csv` or final Merqury metrics, and logs a
  warning when `merqury.sh` exits successfully but no parseable metrics are
  found.
- **Improved** Merqury output organization.  New Merqury runs now use a real
  per-assembly directory and prefix such as `merqury/canu/canu.qv` and
  `merqury/canu/canu.completeness.stats`, while still reading legacy flat
  outputs such as `merqury/canu.qv` from older runs.

---

## [1.3.0] — 2026-04-24

### Overview

Version 1.3.0 is a comprehensive overhaul addressing platform compatibility,
scoring, BUSCO lineage defaults, Merqury integration, and assembly safety.

### Architecture and platform changes

- **Fixed** hifiasm platform compatibility: HiFi-only as primary input.
  ONT ultra-long reads can supplement HiFi via `--ul` but require HiFi as
  primary — not auto-enabled.  CLR disabled.
- **New** assemblers: LJA (La Jolla Assembler, HiFi-only) and raven (all
  platforms) added to default comparison.  MBG available via Bioconda or
  source builds.
- **New** Step 0 — Input QC: validates FASTQ exists, estimates coverage,
  warns for low coverage per-platform (HiFi <25×, ONT <40×, CLR <50×),
  logs compatible assemblers.  Runs in both full and assembly-only modes.
- **Renumbered** pipeline to 15 public steps (0-14):
  - Steps 1-6: original assemblers (Canu, NextDenovo, Peregrine, IPA, Flye, Hifiasm)
  - Steps 7-9: **new** assemblers (LJA, MBG, Raven)
  - Step 10: normalize + QC comparison (BUSCO + Telomere + QUAST + Merqury → assembly_info.csv)
  - Step 11: build telomere pool (pairwise quickmerge + structural validation)
  - Step 12: backbone selection + telomere-aware refinement
  - Step 13: final QC (BUSCO + Telomere + QUAST + Merqury on refined assembly)
  - Step 14: report + cleanup — auto-selects sub-mode:
    14A (full): final comparison report + cleanup into `final_results/`
    14B (assembly-only): assembly-only comparison + cleanup
- **Unified** assembler lists: all downstream code (BUSCO CSV, QUAST CSV,
  Merqury CSV, assembly_info, backbone selection, chimera check) now imports
  `ALL_ASSEMBLERS` from `utils.py` instead of hardcoding names.  Adding a
  new assembler requires only one change in `utils.py`.
- **Assembly-only mode** (`--assembly-only`) runs Steps 0-10, 14
  (assemblers + normalize/QC + assembly-only report).
  Full mode runs Steps 0-14 (adds telomere pool, refinement, final QC,
  report+cleanup).
- **New** taxon-aware BUSCO lineage defaults: `--taxon fungal` → ascomycota_odb10,
  `--taxon plant` → embryophyta_odb10, `--taxon vertebrate` → vertebrata_odb10,
  `--taxon insect` → insecta_odb10, `--taxon other` → requires explicit `--busco`.
  No more fungal-biased default for non-fungal genomes.
- **New** Merqury auto-enabled for ALL platforms: when `merqury.sh` + `meryl`
  are installed, Merqury is enabled by default for HiFi, ONT, and CLR.  TACO
  builds a reads `.meryl` database automatically from input reads, runs Merqury
  on every assembler output (Step 10), and on the final refined assembly
  (Step 13).  For non-HiFi platforms, a warning is logged: QV values may
  underestimate true quality, so Merqury completeness and relative QV ranking
  across assemblers should be interpreted cautiously.  Disable with
  `--no-merqury`.
- **Improved** backbone scoring with taxon-aware weights: BUSCO_S × 1000 -
  BUSCO_D × taxon_penalty + MerquryComp × 200 + MerquryQV × 20 + T2T ×
  taxon_t2t + single × 150 + log10(N50) × taxon_n50 - contigs × taxon_frag -
  size_deviation_penalty.  Insect weights added.  Size deviation penalty
  discourages assemblies far from expected genome size.
- **Improved** replacement safety: candidates must pass sufficient target
  coverage, query coverage, BUSCO validation (C/D/M), size sanity, and
  read coverage support before replacing backbone.
- **New** Step 12 quickmerge structural validation: parent alignment identity,
  query coverage, union coverage, unexplained gap check.  Decision table
  written to `quickmerge_validation.tsv` and `telomere_pool_decisions.tsv`.
- **New** `--merqury-k` flag for configurable Merqury k-mer size (default
  `auto`; use a fixed integer such as 21 or 31 for cross-run comparability).
  Auto-built databases are named `reads.k{K}.meryl`.
- **New** Step 12L "do no harm" comparison: after refinement, compares final
  vs backbone for size, telomere count, and genome size deviation.  If quality
  degraded, saves both assemblies + `refinement_warning.txt`.

### Key changes

- **Fixed** pool T2T upgrade coverage direction bug: the upgrade check now
  verifies BOTH query coverage (pool contig) AND target coverage (backbone
  contig ≥80%).  Previously only checked query coverage, allowing a 1.06M
  pool T2T to "upgrade" a 1.73M backbone contig (losing 667K bp of content
  and ~77 BUSCO genes).
- **New** D-aware duplicate filter (12F2) with three-tier logic for novel
  pool T2T contigs:
  - Overlaps Tier 1 (T2T) backbone: reject (pure duplicate).
  - Overlaps Tier 2 (non-T2T) backbone at ≥80% target coverage: upgrade
    (replace backbone with T2T — better telomere evidence).
  - Overlaps Tier 2 at 50–80% target coverage: **read-coverage diagnostic**
    — maps input reads with a platform-specific minimap2 preset and compares
    median coverage in the T2T-covered region vs the uncovered region.  If uncovered region
    has < 30% of covered region's coverage → backbone is chimeric, replace
    with T2T.  If normal → backbone is real, reject novel (would increase D).
    Configurable via `CHIMERIC_COV_RATIO` (default 0.30).
  - Overlaps Tier 2 at < 50%: reject (insufficient overlap).
  - No overlap: add as genuinely novel (with optional BUSCO D check).
- **New** diagnostic logging for un-upgraded Tier 2 backbone contigs
  (12D2b): reports which Tier 2 contigs have no full T2T upgrade and
  identifies partial T2T hits, flagging potentially chimeric backbone.
- **Disabled** backbone self-dedup by default.  purge_dups handles haplotig
  removal using read-coverage evidence.  Re-enable with `SELFDEDUP_ENABLE=1`.
- **Improved** purge_dups taxon strategy: fungi now use two-round purging
  (`-2`) for aggressive duplicate detection in haploid genomes.  Vertebrate,
  animal, and insect presets use two-round purging for heterozygous
  haplotig/overlap cleanup.  Plants use conservative single-round settings.
  TACO keeps the unpurged assembly when `dups.bed` is empty or when the
  purged result fails taxon-aware size safety checks.  Coverage cutoffs are
  logged.  Override with `PURGE_DUPS_CALCUTS`.
- **Fixed** NextPolish2 v0.2.2: requires sorted BAM as first argument.
  TACO now maps reads with minimap2, sorts with samtools, passes BAM correctly.
- **Fixed** final BUSCO caching: always clears stale results before rerun
  (now part of Step 13 final QC).
- **Fixed** GFF provenance: backbone contigs show clean names (purge_dups
  suffix stripped), `source_type=assembler`; provenance TSV lookup checks
  `final_results/` fallback.
- **Fixed** cleanup file move: uses explicit remove-then-copy to avoid silent
  failures.
- **New** final assembly coverage QC (Step 12K): maps reads to the final
  assembly and computes sliding-window coverage (default 5 kb window).  Detects
  zero-coverage gaps, very low coverage regions (< 15% of global median), and
  mixed low-coverage windows.  Three output files in `final_results/`:
  - `coverage_summary.tsv` — per-contig stats: median, mean, min, max
    coverage, zero-coverage bp, low-coverage bp.
  - `weak_regions.tsv` — every flagged window with coordinates, contig total
    length, source assembler, source type (assembler/quickmerge), window
    median, global median, coverage ratio, and flag (ZERO_GAP, VERY_LOW,
    LOW, MIXED_LOW).
  - `weak_regions.gff3` — GFF3 annotation track for genome browsers (IGV,
    JBrowse).  Each record has `type=coverage_warning`, score = window median
    coverage (for color-coding), and attributes: `flag`, `window_median`,
    `global_median`, `ratio`, `source_assembler`, `source_type`, `description`.
    Load alongside `final.merged.fasta` and `final.merged.provenance.gff3`
    to visualize both provenance and coverage quality in one view.
  Configure window size with `COV_QC_WINDOW` (default 5000) and low-coverage
  threshold with `COV_QC_LOW` (default 5).
- **Fixed** BUSCO trial `busco_available`: no longer requires `--busco` flag.
- **Fixed** `_write_candidates_tsv` KeyError for pool T2T upgrade candidates.
- **Fixed** post-upgrade dedup: protects all backbone contigs, not just
  telomere-bearing ones.

### Environment variables (new)

| Variable | Default | Description |
|----------|---------|-------------|
| `NOVEL_DUP_COV` | 0.80 | Min query coverage to consider novel T2T a duplicate |
| `NOVEL_DUP_ID` | 0.90 | Min identity for duplicate detection |
| `NOVEL_UPGRADE_TCOV` | 0.80 | Min backbone coverage for T2T upgrade |
| `NOVEL_MAX_D_RISE` | taxon default | Max BUSCO D rise allowed for novel additions |
| `PARTIAL_T2T_MIN_TCOV` | taxon default | Min backbone coverage for coverage-guided partial T2T replacement |
| `PARTIAL_T2T_MIN_BP` | taxon default | Min T2T contig length for partial replacement testing |
| `PARTIAL_T2T_MIN_QCOV` | taxon default | Min T2T query coverage for partial replacement testing |
| `PARTIAL_T2T_MIN_IDENT` | taxon default | Min identity for partial replacement testing |
| `STEP12_MAX_ACCEPTED` | taxon default | Max accepted rescue candidates in refinement |
| `STEP12_MIN_BP_RATIO` | 0.90 | Min donor/backbone bp ratio for replacement candidates |
| `STEP12_BUSCO_TRIAL_TIMEOUT` | 43200 | Timeout in seconds for each Step 12 BUSCO trial attempt; 0 disables |
| `STEP12_BUSCO_ALLOW_DOWNLOAD` | 0 | Set to 1 to permit online BUSCO lineage fallback during Step 12/final BUSCO |
| `STEP12_SKIP_BUSCO_TRIAL` | 0 | Set to 1 to use structural rescue validation only |
| `STEP12_MAX_BUSCO_C_DROP` | taxon default | Max BUSCO C% drop allowed during trial validation |
| `STEP12_MAX_BUSCO_M_RISE` | taxon default | Max BUSCO M% rise allowed during trial validation |
| `STEP12_MAX_BUSCO_D_RISE` | taxon default | Max BUSCO D% rise allowed during trial validation |
| `CHIMERIC_COV_RATIO` | 0.30 | If uncovered region coverage < this × covered → chimeric |
| `PURGE_DUPS_CALCUTS` | auto | Override calcuts coverage thresholds |
| `PURGE_DUPS_MODE` | auto | Override purge_dups profile: auto, single, two-round, or skip |
| `PURGE_DUPS_EXTRA_OPTS` | empty | Extra options appended to the purge_dups command |
| `PURGE_DUPS_GET_SEQS_OPTS` | -e | Options for get_seqs; default removes only end duplications |
| `PURGE_DUPS_MAX_BP_DROP` | taxon default | Reject purged output if bp loss exceeds this fraction |
| `PURGE_DUPS_MIN_EXPECTED_RATIO` | taxon default | Reject purged output below this fraction of expected genome size |
| `SELFDEDUP_ENABLE` | 0 | Set to 1 to re-enable backbone self-dedup |
| `COV_QC_WINDOW` | 5000 | Sliding window size (bp) for coverage QC |
| `COV_QC_LOW` | 5 | Reads below this depth counted as low-coverage |
---

## [1.2.0] — 2026-04-15

### Overview

Version 1.2.0 is a major refactor of Step 12 adopting a T2T-first assembly
philosophy.  T2T contigs from all assemblers form the primary foundation;
backbone contigs serve as gap-fill.  Redundans is removed from the pipeline.
Step 12 now performs telomere-aware rescue with donor telomere verification,
BUSCO trial validation, aggressive dedup of non-telomeric backbone contigs,
and self-dedup.  purge_dups replaces Redundans (taxon-aware, ploidy-safe),
and platform-aware polishing (skip for HiFi, Medaka for ONT, Racon for CLR)
is default.  Telomere detection is taxon-aware via `--taxon`, and `--motif`
is an optional override.

---

### Step 12 — T2T-first telomere-aware refinement with BUSCO trial validation

- **Removed** Redundans from the default Step 12 workflow and from `taco-env.yml`.
- **New** T2T-first assembly philosophy: T2T contigs from all assemblers
  form the primary foundation.  Backbone contigs serve as gap-fill only
  for chromosomal regions not covered by T2T contigs.
- **New** backbone telomere classification (12D3): after initial dedup, all
  remaining backbone contigs are classified by telomere status.  Contigs with
  telomere signal are protected from aggressive dedup.
- **New** aggressive non-telomeric dedup (12D4): backbone contigs lacking
  telomere support that overlap the T2T pool at 70%/85% are removed.  Configurable
  via `AGGR_NONTELO_COV` and `AGGR_NONTELO_ID` environment variables.
- **New** non-telomeric self-dedup (12D5): when two non-telomeric backbone contigs
  overlap at 80%/90%, the shorter one is removed.  Telomere-bearing contigs
  are always kept.  Configurable via `SELFDEDUP_COV` and `SELFDEDUP_ID`.
- **New** donor telomere verification: rescue candidates must carry verified
  telomere signal.  Non-telomeric donors are rejected regardless of structural
  alignment quality, enforcing the T2T-first principle.
- **New** structural rescue candidate screening with detailed per-hit metrics:
  identity, aligned bp, backbone/donor coverage, extension, terminal touch,
  length gain, and composite structural score.
- **New** BUSCO trial validation: each plausible candidate is tested by building
  a trial assembly and running BUSCO.  Rejection thresholds are taxon-aware.
- **New** telomere-aware genome-size pruning (12J): telomere-bearing contigs
  are never pruned, only non-telomeric fragments are removed when assembly
  exceeds the genome size budget.
- **New** output files: `single_tel.replaced.debug.tsv` (all hits),
  `single_tel.candidates.tsv` (plausible candidates),
  `rescue_rejection_summary.txt`, `rescue_trial_summary.tsv`.
- **New** provenance GFF3: `final.merged.provenance.gff3` documents each
  contig's full provenance with attributes: `source_assembler`, `role`
  (backbone / upgrade_donor / novel_t2t), `assembler_contig` (original
  assembler contig name before telomere-pool renaming), `source_type`
  (assembler / quickmerge), `replacement_class`, `replaced_contig`, and
  `description` (human-readable provenance summary showing exactly which
  assembler contig replaced which backbone contig and how).
- **New** quickmerge region-level provenance: for contigs derived from
  quickmerge, the GFF3 now includes `qm_assembler1` and `qm_assembler2`
  attributes on the contig-level record, plus child records (type=region,
  linked via Parent) showing which assembler contributed each genomic
  region.  Each region record has `source_assembler` and `assembler_contig`
  tracing to the root original assembler contig.  This enables per-base
  provenance tracing for merged contigs (e.g., "Region 1-500000 from canu
  contig 'tig00001', Region 400000-900000 from flye contig 'contig_3'").
  Region boundaries are determined by minimap2 alignment of each validated
  quickmerge contig against both source assemblers in Step 10.
- **New** `pool_contig_provenance.tsv`: comprehensive provenance map saved
  in Step 10, tracking each pool contig's source assembler, original name,
  and whether it came from an assembler or quickmerge.  For quickmerge
  contigs, extended columns `qm_assembler1`, `qm_assembler2`, and
  `qm_regions` (semicolon-delimited `start-end:assembler:contig` entries)
  record region-level origin.  Used by Step 12 for accurate GFF annotation.
  Moved to `final_results/` during cleanup.
- **Fixed** pool T2T upgrade coverage check: now verifies BOTH query coverage
  (pool contig) AND target coverage (backbone contig ≥80%).  Previously only
  checked query coverage, allowing a small pool contig to "upgrade" a much
  larger backbone contig, losing unique content and BUSCO genes.
- **New** D-aware duplicate filter for novel T2T additions (12F2): before
  adding a "novel" pool T2T contig, aligns it against the current backbone.
  Three outcomes: (a) if it overlaps a Tier 2 (non-T2T) backbone contig, it
  REPLACES the backbone (upgrade — T2T is better than non-T2T); (b) if it
  overlaps a Tier 1 (already T2T) backbone contig, it's rejected as a
  duplicate; (c) if no significant overlap, it's added as genuinely novel.
  Optional BUSCO D check rejects additions that increase D beyond the taxon
  threshold.  Configurable: `NOVEL_DUP_COV` (default 0.80), `NOVEL_DUP_ID`
  (default 0.90), `NOVEL_MAX_D_RISE` (default: taxon D threshold).
- **Disabled** backbone self-dedup by default: even conservative thresholds
  removed backbone contigs with unique BUSCO genes.  purge_dups at 12H
  handles haplotig removal more safely.  Re-enable with `SELFDEDUP_ENABLE=1`.
- **Improved** purge_dups taxon-aware strategy: fungi/haploid genomes now use
  two-round purging (`-2`) with upstream-default chaining lengths; vertebrate,
  animal, and insect genomes use two-round purging for heterozygous duplicate cleanup.
  Plants use conservative single-round settings to preserve homeologs.  TACO
  runs `get_seqs -e` by default, treats an empty `dups.bed` as a valid
  no-duplicates result, and rejects likely over-purged output before replacing
  the assembly.  If the purged assembly moves an overlarge input closer to the
  expected genome size without dropping below the floor, the size drop is
  accepted.  Coverage cutoffs are logged for debugging.  Override with
  `PURGE_DUPS_CALCUTS` env var.
- **Improved** partial strict-T2T handling: substantial taxon-aware partial
  T2T hits now trigger a read-coverage diagnostic on the unmatched backbone
  sequence. Low uncovered coverage promotes the T2T contig to replace a
  duplicated/chimeric backbone contig; normal or inconclusive coverage keeps
  the backbone to preserve BUSCO content.
- **Fixed** NextPolish2 v0.2.2 invocation: requires sorted BAM as first
  argument (`nextPolish2 -t N reads.sorted.bam genome.fa k21.yak k31.yak`).
  TACO now maps HiFi reads with `minimap2 -ax map-hifi`, sorts with
  `samtools sort`, and passes the BAM correctly.  Requires `samtools`.
- **Fixed** final BUSCO caching: always clears stale `busco/final/` results
  before rerunning, preventing silent reuse of outdated metrics.
- **Fixed** cleanup file move: uses explicit remove-then-copy instead of
  `shutil.move` to avoid silent failures when destination already exists.
- Structural screening thresholds: identity >= 0.85, aligned_bp >= 8000,
  cov_backbone >= 0.60, cov_donor >= 0.50, extension >= 1000,
  terminal touch window = 500 bp.  Configurable via environment variables.

### Chimera detection improvements

- **Improved** chimera safety now uses a two-strategy approach:
  1. **Size gate** (existing): contigs > 1.5× the largest individual assembler
     contig are flagged.
  2. **Cross-assembly mapping** (new): each protected contig is aligned against
     all other assembler outputs via minimap2.  A contig not well-covered
     (≥60%) by any other assembler's individual contig is flagged as a
     potential chimera.  Configurable via `CHIMERA_MIN_CROSS_COV` env var.
  Contigs flagged by either strategy are removed from the protected pool.

### Post-dedup BUSCO safety check

- **New** after strict dedup (12D), the pipeline runs BUSCO on the combined
  assembly (protected T2T + remaining backbone) and compares to the backbone
  alone.  If BUSCO C drops > 3% (configurable via `DEDUP_MAX_BUSCO_C_DROP`),
  a prominent warning is logged with remediation suggestions (raise
  `PROTECT_COV`/`PROTECT_ID` thresholds to be more conservative).
- **Improved** `_filter_redundant_to_protected()` now logs each backbone
  contig it removes (name, length, coverage, identity) for easier debugging.

### Post-refinement stack

- **New** purge_dups runs by default after final combine (skip with
  `--no-purge-dups`).  Replaces Redundans for haplotig/overlap cleanup.
  Taxon-aware: uses `-2` (two-round) for vertebrate/animal/plant only;
  single-round for fungal/insect/other to avoid over-purging small genomes.
  Emits a warning for plant genomes (potential polyploid risk).
- **New** automatic polishing runs by default after purge_dups (skip with
  `--no-polish`).  Platform-aware strategy:
  - `--platform pacbio-hifi`: NextPolish2 (k-mer-based polishing via yak
    k-mer databases at k=21 and k=31).  Safe and effective for HiFi data.
    Requires `nextpolish2` + `yak` installed; warns if missing.
  - `--platform nanopore`: Medaka (neural-network polisher) preferred;
    falls back to Racon if Medaka is not installed.
  - `--platform pacbio` (CLR): Racon.

### Taxon-aware telomere detection

- **New** `--taxon` CLI flag: `vertebrate`, `animal`, `plant`, `insect`,
  `fungal`, or `other` (default).  Sets motif-family priors automatically.
- **New** motif families: plant (TTTAGGG) and insect (TTAGG) added to
  `telomere_detect.py`.
- `--motif` is now documented as an optional override, not the recommended
  default.  For fungi and unknown taxa, forcing `--motif` may miss true
  telomeres.

### Assembler platform compatibility

- **Fixed** hifiasm: now correctly skipped for `--platform nanopore` and
  `--platform pacbio` (CLR).  hifiasm only supports HiFi reads as primary
  input.  Previously the pipeline attempted to run hifiasm with an `--ont`
  flag that does not exist.
- **New** taxon-aware backbone scoring weights: fungi penalise BUSCO
  duplicates more heavily and reward T2T contigs; plant/vertebrate reduce
  contig-count penalty for naturally larger assemblies and increase N50 weight.
- **New** taxon-aware BUSCO trial thresholds: fungi 2% C-drop / 0.3%
  M-rise / 2% D-rise (strict); plant 4% / 1.0% / 6% (relaxed for
  polyploidy); vertebrate 3% / 0.5% / 4%.
- **New** BUSCO D-rise (duplicated %) threshold: catches rescue candidates
  that introduce redundant copies of single-copy orthologs.  Configurable
  via `STEP12_MAX_BUSCO_D_RISE` environment variable (`STEP13_*` names remain
  accepted for compatibility with older run scripts).
- **New** two-tier confidence model: Tier 1 (immutable T2T contigs) and
  Tier 2 (editable backbone contigs).  Tier 1 contigs are never replaced
  during rescue unless `--allow-t2t-replace` is explicitly set.
- **New** replacement class tracking: each accepted rescue is assigned one
  of four classes — `fill_missing_end`, `replace_non_telo_backbone`,
  `replace_single_with_better`, or `replace_protected_t2t`.  Classes are
  recorded in `rescue_trial_summary.tsv` and `replaced.ids`.
- **New** taxon-aware rescue limits: fungi 20, vertebrate 10, plant 8,
  other 15.  Prevents runaway replacement in complex genomes.
- **New** taxon-aware 12D4/12D5 dedup thresholds: fungi 70%/85% (12D4)
  and 80%/90% (12D5); plant/vertebrate 85%/92% and 90%/95%; other
  75%/88% and 85%/92%.
- **New** telomere-evidence safety check: `replace_single_with_better`
  candidates are rejected if telomere evidence weakens at either end
  after replacement.
- **New** suspicious size-drop check: candidates causing >30% size drop
  are rejected.

### Version detection improvements

- **Improved** `version.txt` generation: better handling of tools that print
  usage/help to stderr (e.g., seqtk, bwa).  Searches for version-like lines
  in stderr when stdout is empty.
- **New** `purge_dups`, `racon`, `medaka`, and `merqury.sh` added to the
  version.txt tool list.

### Configurability improvements

- **New** `FLYE_ONT_FLAG` environment variable: override Flye read-type flag
  for ONT reads (default `--nano-hq` for Q20+ data; set to `--nano-raw`
  for older pre-Q20 ONT reads).
- **New** `MEDAKA_MODEL` environment variable: override Medaka polishing
  model (default `r1041_e82_400bps_sup_v4.3.0` for R10.4.1 SUP; set to
  `r941_min_sup_g507` for older R9.4.1 data).
- **New** `AGGR_NONTELO_COV`, `AGGR_NONTELO_ID`: thresholds for aggressive
  non-telomeric dedup (default 0.70/0.85).
- **New** `SELFDEDUP_COV`, `SELFDEDUP_ID`: thresholds for non-telomeric
  self-dedup (default 0.80/0.90).
- **New** `CHIMERA_MIN_CROSS_COV`: minimum cross-assembly coverage for
  chimera safety mapping check (default 0.60).
- **New** `DEDUP_MAX_BUSCO_C_DROP`: maximum tolerated BUSCO C drop after
  dedup before a warning is issued (default 3.0%).

### CLI changes

- Added `--taxon` parameter.
- Added `--no-purge-dups` to skip purge_dups.
- Added `--no-polish` to skip polishing.
- Added `--allow-t2t-replace` to permit rescue donors to replace
  immutable Tier 1 (T2T) contigs.  Disabled by default for safety.
- Updated `--motif` help text to indicate it is an override.
- Updated `--reference` help text (Redundans reference removed).

### Environment changes

- Removed `redundans` from `taco-env.yml`.
- Added `purge_dups`, `racon`, and `medaka` to `taco-env.yml`.
- Added `nextpolish2` and `yak` to `taco-env.yml` for default HiFi polishing.

### Taxon-aware score windows

- Default telomere score window is now set by taxon: 300 bp for fungi
  (short telomere arrays), 1000 bp for plant/vertebrate (longer arrays),
  500 bp for others.

### Version bumps

- `__init__.py`: 1.1.0 → 1.2.0
- `setup.py`: 1.1.0 → 1.2.0
- `pipeline.py`: TACO-1.0.0 → TACO-1.2.0
- `cli.py`: version string updated to v1.2.0

---

## [1.1.0] — 2026-04-14

### Overview

Version 1.1.0 addresses a systematic BUSCO duplication problem (D ≈ 57.5% in
the final merged assembly) introduced by chimeric contigs from quickmerge
entering the protected T2T pool.  The release introduces telomere-aware
validated quickmerge in Step 10, assembler quality weighting in clustering,
and corrects the redundancy-filter thresholds in Step 12 that had been
inadvertently lowered.

---

### Step 10 — Telomere-aware validated quickmerge (`steps.py`)

**Problem (root cause of high BUSCO D):** The original `allmerged_telo.fasta`
pool was built by concatenating all pairwise quickmerge outputs together with
the original assembler `.telo.fasta` files.  Quickmerge can join two contigs
from different chromosomes into a single chimeric contig; if that chimera
happens to acquire telomere signal on both ends (e.g. by being flanked by
repeat-rich ends from either input), it classifies as `strict_t2t` and enters
the protected pool.  When the chimera then displaces genuine backbone contigs,
both chromosomes' gene sets appear duplicated — inflating BUSCO D.

**Fix — telomere-aware validated merge:** A new `_validate_quickmerge_t2t()`
function runs after each pairwise quickmerge call and accepts a merged contig
into the pool only when both criteria are met:

1. **Telomere proof.** The merged contig must score as `strict_t2t` (telomere
   signal on both ends, both raw scores ≥ 0.25).  This is positive evidence
   that the join actually rescued a missing telomere end rather than simply
   concatenating two unrelated sequences.
2. **Length sanity.** Merged contig length ≤ 1.3× max(input contig lengths).
   A legitimate join of two overlapping single-end contigs that represent the
   same chromosome (the match–mismatch–match gap-resolution pattern) produces
   a contig close in length to the longer input.  A chimeric join of two
   different chromosomes produces a contig ≈ sum(input lengths), i.e. roughly
   2× — this is rejected.

The pool is thus: original assembler `.telo.fasta` files + validated
quickmerge T2T contigs (tagged `_qm_validated`).  Quickmerge's purpose of
joining a left-telomere contig from one assembler with a right-telomere contig
from another is fully preserved; only the chimeras are excluded.

---

### Step 10 — Assembler BUSCO D quality weighting for clustering (`steps.py`)

**Problem:** Clustering representative selection used telomere score alone as
the tiebreaker between contigs in the same chromosome group.  A contig from a
high-duplication assembler (e.g. canu BUSCO D = 74.9%) could win over a contig
from a clean assembler (e.g. peregrine BUSCO D = 10%) if it had marginally
higher telomere scores, carrying duplicated gene content into the protected
pool.

**Fix:** Before calling `cluster_and_select`, the code now reads
`assemblies/assembly_info.csv` to build an assembler → BUSCO D map.  Each
contig's effective clustering score is multiplied by a quality weight derived
from its source assembler's BUSCO D:

```
quality_weight = max(0.25,  1.0 − busco_d / 133.0)
```

A contig from peregrine (D = 10%) receives weight ≈ 0.92; one from canu
(D = 74.9%) receives weight ≈ 0.44.  The clustering therefore strongly prefers
representatives from lower-duplication assemblers when telomere scores are
similar.  The assembler source is inferred by matching pool contig sequence
lengths back to the original `.telo.fasta` files.

---

### Step 10 — T2T contigs no longer stripped by `_fasta_clean_contained` (`steps.py`)

**Problem:** `_fasta_clean_contained` was being called on the T2T pool before
clustering.  Because T2T contigs on different chromosomes can share repetitive
subtelomeric elements, the 30% coverage threshold occasionally matched a T2T
contig against a longer contig from a different chromosome and removed it.  In
one observed run, 16 T2T contigs entering Step 10 became 15 in the protected
pool (1 chromosome lost its T2T representative).

**Fix:** `_fasta_clean_contained` is no longer called on `t2t.fasta`.  The
minimap2 clustering (`cluster_and_select` with 95% identity / 85% query
coverage) already deduplicates genuine same-chromosome copies; the clean step
added no benefit and caused cross-chromosome false positives.

---

### Step 10 — Per-step log files (`pipeline.py`)

**Problem:** All step output was written only to the terminal; there was no
persistent per-step log for post-run debugging.

**Fix:** Added `TeeWriter` class in `pipeline.py`.  The step execution loop
now redirects both stdout and stderr through `TeeWriter` so that every step
writes a file at `logs/step_N.log` in addition to the console.

---

### Step 9 — Telomere metric name consistency across steps (`steps.py`)

**Problem:** Step 9 wrote rows named `"Telomere double-end contigs"`,
`"Telomere single-end contigs"`, and `"Telomere supported contigs"` to
`assembly_info.csv`.  Steps 12, 14, and 16 looked up `"Telomere strict T2T
contigs"` and `"Telomere single-end strong contigs"`.  The name mismatch
caused all telomere metrics to read as 0 in backbone scoring and final report
generation.

**Fix:** All steps now use a single consistent set of names:

| Tier | Row label in `assembly_info.csv` |
|---|---|
| Both ends ≥ 0.25 | `Telomere strict T2T contigs` |
| One end ≥ 0.25 | `Telomere single-end strong contigs` |
| Any end ≥ 0.08 | `Telomere-supported contigs` |

---

### Step 9 — Telomere detection for all assemblers, not only hifiasm (`steps.py`)

**Problem:** Hybrid telomere detection in Step 9 was running the full
`detect_telomeres()` pipeline only on the hifiasm output.  All other
assemblers fell through to a placeholder path and produced empty
`.telo.fasta` files, making their telomere metrics 0 in all downstream steps.

**Fix:** The detection loop now iterates over every assembler result FASTA
present in `assemblies/` and calls `detect_telomeres()` for each.

---

### Step 12 — Redundancy-filter thresholds reverted to 95 % / 95 % (`steps.py`)

**Problem:** In a previous edit, the `_filter_redundant_to_protected` call
thresholds in Steps 12D and 12F were lowered from 95%/95% to 60%/90%.  The
intent was to remove backbone contigs that partially overlapped a T2T contig,
but the effect was the opposite: more clean backbone contigs (e.g. from
peregrine, D = 10%) were dropped and replaced by T2T pool contigs that came
from high-D assemblers, driving BUSCO D up rather than down.

**Fix:** Both dedup passes (12D backbone-vs-protected and 12F
rescued-vs-protected) now use `cov_thr = 0.95, id_thr = 0.95`, matching the
original TACO.sh default.  These conservative thresholds only drop backbone
contigs that are near-identical to a protected T2T contig; genuine backbone
contigs from independent chromosomes are kept.  The thresholds remain
overridable via `PROTECT_COV` and `PROTECT_ID` environment variables.

---

### Step 12 — BUSCO D penalty added to backbone scoring (`steps.py`)

**Problem:** The smart backbone scoring formula rewarded BUSCO single-copy
completeness but did not penalise duplication.  An assembler with high BUSCO S
but also high BUSCO D could be selected as backbone, carrying duplicated
content into the final assembly.

**Fix:** BUSCO D is now read from `assembly_info.csv` alongside BUSCO S, and
a penalty term is subtracted:

```
score = BUSCO_S×1000 + T2T×300 + single×150
      + MerquryComp×200 + MerquryQV×20
      − contigs×30 + log₁₀(N50)×150
      − BUSCO_D×500
```

The BUSCO D value and the penalty contribution are written to
`assemblies/selection_debug.tsv` and `assemblies/selection_decision.txt` for
transparency.

---

### Step 12 — Chimera safety check on protected pool (`steps.py`)

**Problem:** Even with the validated-merge approach in Step 10, an additional
guard is warranted against any abnormally long contigs that may have entered
the protected pool.  The earlier guard used 1.5× the median backbone contig
length, which was too aggressive for organisms where chromosome sizes vary
significantly (legitimate longer chromosomes would be incorrectly flagged).

**Fix:** The chimera safety check now uses 1.5× the **maximum individual
contig length across all assembler `.telo.fasta` files** as the threshold.
Because this reference length comes from actual chromosome-scale contigs (the
largest chromosomes observed in any assembler output), 1.5× catches genuine
2× chimeras while leaving all legitimate T2T contigs — including those from
larger chromosomes — untouched.

---

### Step 12 — Redundans integration for reduction, scaffolding, and gap closing (`steps.py`, `taco-env.yml`)

**Problem:** The custom minimap2-based fragment removal (Pass 2 at 50%/90%)
only catches backbone contigs that partially align to a protected T2T contig.
It misses redundant heterozygous pairs among the backbone contigs themselves,
and it cannot join fragmented backbone contigs or close assembly gaps.
Additionally, `redundans.py` was listed in `TACO.sh`'s requirements check
(line 772) but was never actually called and was missing from `taco-env.yml`.

**Fix — Redundans on the full combined assembly (Step 12G2):**

Redundans runs after the final combine (protected T2T + backbone) so it can
see both the T2T chromosomes and the surviving backbone fragments.  An earlier
design ran Redundans reduction on the backbone alone (before the combine), but
this produced "Nothing reduced!" because the backbone fragments are redundant
*against the T2T contigs*, not among themselves.

Step 12 now proceeds:

1. **12D Pass 1** — strict minimap2 dedup (95%/95%) removes near-identical
   backbone contigs.
2. **12D Pass 2** — minimap2 fragment removal (50%/90%) removes backbone
   fragments that partially overlap T2T chromosomes.
3. **12E–12G** — telomere rescue, post-rescue dedup, final combine.
4. **12G2 — Redundans on the full combined assembly.**  All three Redundans
   stages run in order:
   - **Reduction** — detects and removes heterozygous/duplicate contigs across
     the entire assembly (identity ≥ 0.51, overlap ≥ 0.80; configurable via
     `RED_IDENTITY` / `RED_OVERLAP` env vars).
   - **Scaffolding** — joins fragments using the original long reads (`-l`).
   - **Gap closing** — fills gaps created during scaffolding.
   The minimap2 preset is auto-selected from `--platform` (map-hifi, map-ont,
   map-pb).
5. **12H** — genome-size-aware pruning (safety net).

**Fallback:** If `redundans.py` is not installed, step 12D keeps the minimap2
fragment removal, and Redundans scaffolding/gap-closing is skipped.

**Reference-guided vs de novo mode:**  When `--reference` / `-ref` is
provided, it is passed to Redundans as `-r` for reference-guided reduction
and scaffolding.  For pure de novo runs (no `--reference`), the reference
is skipped entirely — Redundans uses only HiFi/long reads.

**CLI rename:**  `--fasta` renamed to `--reference` / `-ref` to clarify its
role as a reference genome (not an "external" assembly).  The assembler name
`external` is now `reference` in all CSV outputs, file paths, and comparison
tables (`assembly_info.csv`, `selection_debug.tsv`, `final_result.csv`, etc.).
The internal attribute is `runner.reference_fasta` (was `external_fasta`).
The README documents this dual-mode behaviour.

`taco-env.yml` now includes `redundans` from bioconda.

---

### Clustering — query-only coverage and telomere score tiebreaker (`clustering.py`)

**Problem:** `parse_paf_and_cluster` computed coverage as
`min(query_cov, target_cov)`.  When a shorter contig was fully contained
within a longer one, the target coverage was small and the pair was not
clustered together, leaving genuine duplicates as separate representatives.

**Fix:** Coverage is now computed on the query contig only
(`(qe − qs) / qlen`), consistent with TACO.sh.  Additionally:
- `seq_names` parameter added so singletons (contigs with no alignments)
  always appear as their own single-member cluster.
- `scores` parameter added to `cluster_and_select`; when provided, the
  cluster representative is chosen by `(telomere_score DESC, length DESC,
  name ASC)` rather than length alone.  This is used in Step 10 with the
  BUSCO-weighted scores described above.

---

## [1.0.0] — 2026-04-08

### Overview

Version 1.0.0 is the first stable public release of TACO.  The entire
pipeline has been converted from a monolithic 2,620-line Bash script
(`TACO.sh`) into a proper, installable Python package (`taco/`).
The scientific logic is identical to v0.5.6 — what changed is how that
logic is packaged, invoked, and maintained.

### Summary of changes

- **MAJOR** — Entire pipeline rewritten from Bash to a proper Python package
  (`taco/`); `pip install -e .` registers the `taco` console-script entry point.
- **MAJOR** — Dependency `funannotate` removed; `funannotate sort` replaced by
  `_fasta_sort_minlen()` and `funannotate clean` by `_fasta_clean_contained()`,
  both pure Python + minimap2.
- **FIX** — Bash function-ordering bug (`check_single_env_requirements` called
  before it was defined) eliminated by the Python conversion.
- **FIX** — All embedded Python heredocs extracted into proper module functions;
  no more `python3 <<'EOF' … EOF` patterns.
- **FIX** — All `awk`/`sed` pipelines replaced with Python `csv`, `re`, and
  string methods — no GNU vs BSD `awk` discrepancies.
- **CHANGE** — Assembler steps 1–6 are now non-fatal: binary absence or
  non-zero return code logs `[warn]` and continues to the next assembler.
- **CHANGE** — Platform-specific assembler flags (`-pacbio-hifi`, `--nano-hq`,
  etc.) now selected from `ASSEMBLER_PLATFORMS` dict based on `--platform`.
- **FIX** (canu) — `openjdk>=11` added as explicit conda dependency to prevent
  `undefined symbol: JLI_StringDup` from bioconda dev builds.
- **FIX** — GitHub repository URL corrected from `ysun-fieldmuseum/TACO` to
  `yksun/TACO` in README shields, setup.py, and all internal references.

---

### Architecture — Bash → Python conversion

#### Why the conversion was done

The original `TACO.sh` worked well as a development script but had
structural problems that made it increasingly fragile:

- **Function-ordering bug.** In Bash, a function must be defined before it
  is called at the top level.  `TACO.sh` called
  `check_single_env_requirements` at line 623 but defined it at line 798,
  causing `command not found` every time the script was sourced on a fresh
  environment.  Python has no such restriction — functions and classes can
  be referenced freely anywhere in a module.

- **Embedded Python heredocs.** Many steps contained Python programs
  written as Bash here-documents (`` python3 <<'EOF' … EOF ``).  These were
  difficult to edit, impossible to test in isolation, and produced cryptic
  indentation errors when the surrounding Bash changed.

- **Fragile awk/sed pipelines.** Column parsing, CSV construction, telomere
  classification, and BUSCO result aggregation all relied on long
  `awk '{…}'` chains piped through `sed` and `tr`.  These broke silently on
  different `awk` builds (GNU vs BSD), on files with Windows-style `\r\n`
  line endings, and when upstream tool output formats changed slightly.

- **No installability.** Users had to call `bash ~/opt/TACO/TACO.sh -flags`
  and manually manage `PATH` and `PYTHONPATH`.  There was no standard entry
  point and no way to use `taco --help` directly.

#### What the conversion looks like

| Bash pattern | Python replacement |
|---|---|
| `function foo() { … }` at arbitrary positions | `def foo(runner):` in any module; no ordering constraint |
| Top-level variable assignments (`THREADS=30`) | `PipelineRunner` instance attributes (`self.threads`) |
| `export PYTHONPATH=…; python3 <<'EOF' … EOF` | Direct call to a Python function in `taco/steps.py` |
| `awk '{print $3}' file \| sort \| uniq` | Python `csv.reader`, `dict`, and `set` operations |
| `sed 's/\r//'` for Windows line endings | Python `line.rstrip('\r\n')` in every file reader |
| `if [ -f "$f" ] && [ -s "$f" ]` | `os.path.isfile(f) and os.path.getsize(f) > 0` |
| `bash_array=("a" "b"); for x in "${bash_array[@]}"` | Python list: `for x in ["a", "b"]` |
| `assembler_result=$(cat file \| grep …)` | Python string methods and `re` module |
| Assembler step exits with `exit $?` on failure | `check=False` + `result.returncode != 0` → `log_warn` + `return` |
| `funannotate sort -i in -b contig -o out --minlen 500` | `_fasta_sort_minlen(in, out, prefix="contig", minlen=500)` |
| `funannotate clean -i in -p 30 -o out --exhaustive` | `_fasta_clean_contained(in, out, pct_cov=30, exhaustive=True)` |

#### New package layout

```
taco/                    (2,934 lines total)
├── __init__.py          version string
├── __main__.py          entry point for `python3 -m taco` and `taco`
├── cli.py               argparse — all flags for all public steps
├── pipeline.py          PipelineRunner class — logging, run_cmd, version checks
├── steps.py             step functions + FASTA helpers
├── utils.py             ASSEMBLER_PLATFORMS dict, FASTA I/O utilities
├── telomere_detect.py   hybrid telomere detection (MOTIF_FAMILIES, scoring)
├── telomere_pool.py     three-tier pool classification
├── clustering.py        UnionFind, PAF-based contig clustering
├── backbone.py          smart/N50 backbone scoring and selection
└── reporting.py         final CSV report generation

setup.py                 pip install -e . → registers `taco` console script
taco-env.yml             conda environment with all dependencies
run_taco                 optional wrapper for use without pip install
```

#### How to install after the conversion

```bash
conda env create -f taco-env.yml
conda activate taco
pip install -e .          # registers the `taco` command in PATH
taco --help
```

Previously: `bash ~/opt/TACO/TACO.sh -t 30 --fastq reads.fastq -g 12m`
Now:        `taco -t 30 --fastq reads.fastq -g 12m`

---

### Assembler robustness — steps 1–6 now non-fatal

**Previously:** any assembler failure caused the entire pipeline to exit
immediately via `sys.exit(returncode)`.  A broken canu Java runtime would
abort everything before flye or hifiasm even started.

**Now:** each of the six assembler steps (canu, nextDenovo, peregrine, IPA,
flye, hifiasm) follows the pattern:

1. Check that the binary exists (`shutil.which`).  If absent, log a `[warn]`
   with the install command and `return` — the next assembler continues.
2. Run the assembler with `check=False` (no automatic exit on failure).
3. If `returncode != 0`, log `[warn] Step N failed. Skipping.` and `return`.

The pipeline always proceeds to the next step regardless of how many
assemblers fail.

### Platform-specific assembler flags

All assembler steps previously hardcoded `-pacbio-hifi` (or the HiFi
equivalent).  They now select the correct flag from `ASSEMBLER_PLATFORMS`
in `utils.py` based on `--platform`:

| Assembler | pacbio-hifi | nanopore | pacbio (CLR) |
|---|---|---|---|
| canu | `-pacbio-hifi` | `-nanopore` | `-pacbio` |
| nextDenovo | `hifi` (cfg) | `ont` (cfg) | `clr` (cfg) |
| peregrine | ✓ supported | ✗ skipped | ✓ supported |
| IPA | ✓ supported | ✗ skipped | ✗ skipped |
| flye | `--pacbio-hifi` | `--nano-hq` | `--pacbio-raw` |
| hifiasm | default | `--ont` | default |

### Canu — Java runtime fix

The bioconda `canu` package frequently installs a development build
(e.g. `r10515 master +27 changes`) whose bundled Java produces
`undefined symbol: JLI_StringDup` at runtime.

**Fix:** `openjdk>=11` is now an explicit dependency in `taco-env.yml`.
Conda resolves a working JRE before installing canu, preventing the
dynamic-linker failure.  Additionally, the pipeline detects dev builds by
inspecting the version string at runtime and emits a `[warn]` directing
the user to the stable release page before attempting assembly.

### Dependency removed: funannotate

`funannotate` was a ~200 MB annotation suite required only for two FASTA
utility subcommands.  Both are now implemented directly in `steps.py` with
no new dependencies:

**`funannotate sort`** → `_fasta_sort_minlen(infa, outfa, prefix, minlen)`:
- Reads the input FASTA, discards sequences shorter than `minlen` (500 bp),
  sorts remaining sequences by length (longest first), and renames them
  `prefix_1`, `prefix_2`, …  Pure Python, no subprocess.

**`funannotate clean`** → `_fasta_clean_contained(infa, outfa, pct_cov, exhaustive)`:
- Self-aligns with `minimap2 -x asm5 -DP` (minimap2 is already a required
  dependency).  Parses the PAF output: for each alignment where the shorter
  query contig is covered ≥ `pct_cov`% (30%) by a longer target contig, the
  query is marked for removal.  If `exhaustive=True`, repeats until no
  further contigs are removed (up to 20 rounds).  The algorithm is
  functionally identical to `funannotate clean --exhaustive`.

### GitHub repository URL corrected

All URLs, badge links, and internal references now point to the correct
repository `yksun/TACO` (previously `ysun-fieldmuseum/TACO`).  Affected
files: `README.md` (4 shields.io badges + clone URL), `setup.py` (`url=`).

---

## [0.5.6] — 2026-03-xx

- **MAJOR** (Steps 9, 10, 14): Replaced `seqtk telo` exact-motif detection
  with a hybrid scoring system: de novo k-mer discovery + built-in
  MOTIF_FAMILIES (canonical TTAGGG, budding yeast TG1-3, Candida 23-bp) +
  per-end composite scoring.  Three-tier classification:
  strict_t2t (both ends ≥ 0.25), single_tel_strong (one end ≥ 0.25),
  telomere_supported (one end ≥ 0.08).
- **FIX** (Step 12): Backbone scoring reads updated telomere column names
  (`Telomere strict T2T contigs`, `Telomere single-end strong contigs`)
  instead of the old names used in v0.5.x.
- **FIX** (Step 10): Fixed `awk` variable scope bug where `LIST` was passed
  as a positional argument rather than `-v LIST=…`, making it unavailable in
  the `BEGIN` block and producing empty telomere pool extraction.
- **FIX** (Step 16): Final report now uses score-based telomere
  classification from `.telomere_end_scores.tsv` instead of
  position-based heuristics.
- **CHANGE** (Step 12): Backbone scoring switched from BUSCO C% to BUSCO S%
  to avoid rewarding duplicated assemblies.  Updated formula:
  `BUSCO_S×1000 + T2T×300 + single×150 + MerquryComp×200 + MerquryQV×20 − contigs×30 + log₁₀(N50)×150`
- **FIX** (Step 10): Pool classification `.list` files now written directly
  by Python, removing `awk` TSV extraction as a failure point.
- **MAJOR** (Step 10): Strict meaning of `t2t.fasta` preserved — contains
  only true double-end T2T contigs.  `single_tel_best.fasta` holds best
  single-end representatives; `telomere_supported_best.fasta` holds the
  combined optimized pool.
- **MAJOR** (Step 10): Redundancy reduction for single-end telomeric contigs
  via all-vs-all minimap2 clustering and longest-representative selection.
- **MAJOR** (Step 10): Protected telomere-pool priority: strict T2T >
  best single-end > optimized telomere-supported.
- **CHANGE** (Step 18): Added `--assembly-only` mode for benchmarking
  without telomere-aware refinement.
- **DESIGN**: Pipeline now follows telomere-pool optimization + backbone
  refinement rather than repeated structural merging.

---

## [0.5.5] — 2026-03

Identical science to 0.5.6; minor internal fixes to pool classification
and column naming.  See 0.5.6 for full description of the hybrid telomere
detection system introduced in this development series.

---

## [0.5.4] — 2026-02

- **MAJOR** (Step 10): Reworked telomere-pool construction into three
  biologically distinct classes (strict T2T, best single-end, optimized
  telomere-supported).
- **MAJOR** (Step 12): Backbone refinement logic replaces old
  `merge_wrapper`-based final merge; the selected assembler output is used
  as the backbone with telomere-pool replacement.
- **CHANGE** (Step 12): Smart/N50 automatic backbone selection modes added.
- **CHANGE** (Step 12): Optional Merqury integration for assembler ranking.
- **FIX** (Step 12): Relaxed fungal single-end telomere rescue thresholds
  (terminal overhang, alignment identity, coverage).

---

## [0.5.3] — 2026-02

- **FIX** (Step 16): Telomere rescue pool now prioritizes
  `single_tel_best_clean.fasta` before broader fallback sets.
- **CHANGE** (Step 14): Final telomere analysis reports strict T2T,
  single-end, total telomere-supported, protected mode, and rescue counts.
- **CHANGE** (Step 16): Expanded `final_result.csv` to include Merqury
  metrics, rescue counts, selection score, selected assembler, and
  auto-selection mode.

---

## [0.5.2] — 2026-01

- Reworked telomere-pool construction to prioritize strict T2T and select
  best single-end representatives using minimap2 clustering.
- Added `telomere_cluster_summary.tsv` and updated telomere support
  summary for representative-contig selection.
- Updated backbone refinement to use the optimized telomere pool.
- Increased telomere contribution in smart scoring.

---

## [0.5.1] — 2025-12

- **FIX** (Step 12): Unified `--choose` / `--auto-mode` argument parsing
  to prevent auto-mode from being ignored.
- **CHANGE** (Step 12): Merqury pre-selection evaluation added for all
  available assembler outputs when a `.meryl` database is present.
- Added `selection_debug.tsv` and `selection_decision.txt` for transparent
  scoring and reproducibility.
- **FIX** (Step 12): Post-rescue deduplication against protected telomere
  contigs reduces redundancy before final assembly.
- **CHANGE** (Step 16): Final Merqury evaluation for
  `assemblies/final.merged.fasta`; output renamed to
  `final_results/final_result.csv`.

---

## [0.5.0] — 2025-12

- **MAJOR** (Step 10): Three-tier telomere classification from `seqtk telo`
  coordinates: `t2t.fasta` (strict double-end), `single_tel.fasta`
  (single-end), `telomere_supported.fasta` (union).
- **MAJOR** (Step 12): `merge_wrapper`-based final merge replaced by
  backbone-refinement strategy using `protected_telomere_contigs.fasta`.
- Added `protected_telomere_mode.txt` for downstream rescue logic.
- **DESIGN**: Protection/replacement strategy replaces repeated structural
  merging; telomere-supported contigs are preserved first, best assembly
  used as nonredundant backbone.

---

## [0.30] — 2025-11

- Integrated `t2t_list.sh` as an internal function; no external helpers
  needed for Steps 9 and 14.
- **FIX** (Step 12): CLI args passed to Python heredocs correctly
  (no post-heredoc args) — fixes `IndexError` and `Permission denied`.
- Robust project-name parsing for `.fastq.gz` inputs.
- Exposed `PROTECT_COV`/`PROTECT_ID` environment variables for
  redundancy filter (default 0.95/0.95).
- **FIX** (Steps 9 & 14): T2T contigs counted when both ends have
  telomeric signal (start == 0 and end == sequence length).
- **FIX** (Step 12): `awk` CR-stripping made portable
  (`gsub("\\r","")` instead of regex literal).
- **FIX** (Step 7): Removed stray Bash inside Python heredoc.

---

## [0.2.7] — 2025-11-14

- **CHANGE** (Step 12): Normalize `others.fa` headers; run redundans with
  minimap2-based reduction; log contig-reduction statistics.

---

## [0.2.6.8] — 2025-11-12

- **CHANGE** (Step 12): Auto-select assembler with highest N50 from
  `assemblies/assembly_info.csv` when `--choose` is not provided.

---

## [0.2.6.7] — 2025-11-12

- **FIX** (Step 12): Guard against missing/empty `t2t_clean.fasta` —
  skip `merge_wrapper` and pass through the chosen assembler result
  directly when no T2T contigs exist.

---

## [0.2.6.6] — 2025-11-12

- **FIX** (Step 10): Proceed with available FASTAs when some assemblers
  are absent.  For a single available assembly, skip pairwise merge and
  copy directly to `allmerged_telo.fasta`.
