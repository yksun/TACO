# TACO Changelog

All notable changes to TACO are documented here.
Versions follow [Semantic Versioning](https://semver.org/).

---

## [1.0.0] — 2026-04-08

### Overview

Version 1.0.0 is the first stable public release of TACO.  The entire
pipeline has been converted from a monolithic 2,620-line Bash script
(`TACO.sh`) into a proper, installable Python package (`taco/`).
The scientific logic is identical to v0.5.6 — what changed is how that
logic is packaged, invoked, and maintained.

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
