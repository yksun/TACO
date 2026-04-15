# TACO Changelog

All notable changes to TACO are documented here.
Versions follow [Semantic Versioning](https://semver.org/).

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

**Fix — three-stage Redundans integration:**

1. **Redundans reduction (Step 12D, Pass 2).** After the strict 95%/95%
   minimap2 dedup, the surviving backbone contigs are passed to
   `redundans.py --noscaffolding --nogapclosing --minimap2reduce`.  Redundans
   detects and removes heterozygous/duplicate contigs using its internal
   alignment + identity/overlap logic (defaults: identity ≥ 0.51,
   overlap ≥ 0.80).  Thresholds are overridable via `RED_IDENTITY` and
   `RED_OVERLAP` environment variables.

2. **Redundans scaffolding + gap closing (Step 12G2, after final combine).**
   The full combined assembly (protected T2T + backbone) is passed to
   `redundans.py --noreduction` with the original long reads (`-l`).  This
   step joins backbone fragments using read evidence (scaffolding) and fills
   the gaps using the same reads (gap closing).  The minimap2 preset is
   auto-selected from `--platform` (map-hifi, map-ont, map-pb).

3. **Fallback.** If `redundans.py` is not installed, step 12D falls back to
   the minimap2-based fragment removal (50%/90%), and scaffolding/gap-closing
   is skipped with an informational message.

**Reference-guided vs de novo mode:**  When `--fasta` (external) is provided,
it is passed to Redundans as `-r` (reference) for both reduction and
scaffolding.  This enables reference-guided scaffolding where Redundans uses
the reference chromosome structure to order and orient contigs.  For pure de
novo runs (no `--fasta`), Redundans operates without a reference and relies
solely on contig-vs-contig alignment (reduction) and long-read evidence
(scaffolding).  The README documents this dual-mode behaviour.

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
├── cli.py               argparse — all flags for all 18 steps
├── pipeline.py          PipelineRunner class — logging, run_cmd, version checks
├── steps.py             all 18 step functions + FASTA helpers
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
- Updated backbone refinement to use the optimized Step 10 pool.
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
