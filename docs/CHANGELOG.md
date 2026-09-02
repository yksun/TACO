# TACO Changelog

All notable changes to TACO are documented here.
Versions follow [Semantic Versioning](https://semver.org/).

---

## [1.5.2] — 2026-09-02

Five corrections found while re-running v1.5.1 on *Fusarium tricinctum* and
re-reading its logs. None changes the delivered assembly on that data set. Three
change what TACO reports about it, and two stop an in-place rerun from describing
the previous run's files as this run's.

### A metric that some candidates lack no longer penalises them

Merqury does not always produce a QV for every assembly. `score_assembly` treated
a missing value as zero, and the docstring claimed that made the term neutral. It
does not: zero is the *minimum*, so a missing value is the maximum penalty. v1.5.0
made this six times worse by raising `w_merqury_qv` from 20 to 120.

On the real *Fusarium tricinctum* run two of nine candidates had no Merqury QV,
which cost them **6,600 points** against a peer at QV 55 and pushed them from
ranks 2–3 to 6–7.

A metric present for some candidates but not all is not comparable across the
set, so `partially_available_metrics` now detects that case and the term is
dropped for **every** candidate, with a warning naming the metric. Dropping it
for the absentees alone would distort the comparison just as badly in the other
direction. The selected backbone is unchanged on that data set either way, which
is the reassuring part.

The decision record now matches the score. `selection_decision.txt` had listed
`MerquryQV*120` in `score_formula` while the term was excluded from every
candidate's score — exactly the case on the *Fusarium* run, where Flye and
NextDenovo had no QV. `_describe_score_formula` now omits a dropped term and
appends `merqury_qv dropped: missing for some candidates, so excluded for all`
to the formula's notes.

### Cached Merqury output must post-date the assembly it describes

Re-running in place leaves the previous run's `merqury/<label>/` in position, and
the reuse check asked only whether the files existed. Observed on a real rerun: a
120 Mb assembly's QV 76.77 and completeness 88.73 were reported verbatim for the
145 Mb assembly that replaced it. BUSCO already validates the lineage recorded in
its output; Merqury had no equivalent.

- **New `_merqury_cache_is_current`.** Cached `.qv` and `.completeness.stats`
  are reused only when both are newer than the assembly FASTA. Otherwise Merqury
  is recomputed, and the log names the stale file.

### "No unsupported k-mers" is a result, not a missing value

Merqury writes `+inf` in the QV column when every assembly k-mer occurs in the
reads. The parser's sanity guard (1–1000) rejected it and TACO reported `NA`, so
the best possible QV looked like a failed measurement. Observed on the purified
*Fusarium* assembly: 0 error k-mers, against the published reference's 110.

- `_parse_merqury_qv` now returns `inf (no unsupported k-mers)` for `+inf`,
  `inf`, and `Infinity`. The value is written to `final_result.csv` as is, and
  selection treats it as unscoreable, which is correct: it is not comparable to
  a finite QV, and the delivered assembly is not a selection candidate.

### 12L compares against this run's assembly only

The delivered assembly is emitted at 13C and measured at 13D, after 12L runs, so
any `busco/final` present at 12L belongs to an earlier run. On a real rerun 12L
compared this run's backbone against last run's delivered assembly and reported a
plausible-looking BUSCO delta. 12L now uses `busco/final` only when it post-dates
`final.merged.fasta`; otherwise it states that the comparison is deferred to the
final report.

### A no-op purge is described in words

`purge_dups accepted: removed -0.0 points of duplication for -0.0 points of
complete genes` was the log line on a haploid genome with nothing to purge. The
gate now reports signed deltas, and when both round to zero it says
`nothing to purge`.

---

## [1.5.1] — 2026-08-24

### The merge was removing complete genes, and every gate was watching something else

v1.5.0 added a gene-content gate to purge_dups. It was the wrong step. On a real
*Puccinia triticina* run the delivered assemblies still lost gene content in both
representations — BUSCO complete fell 95.6% → 92.5% in primary and **96.7% →
90.4% in full**, which skips purge_dups entirely — so the gate could not have
been responsible.

Tracing it: of the 112 complete BUSCOs lost in full mode, **104 were `Missing`
and only 8 `Fragmented`**. The sequence was gone, not broken at a merge junction.
Missing rose from 2.6% to 8.5%.

The cause is step 12's merge. It rebuilds the assembly from the telomere pool and
drops backbone sequence the pool appears to cover, and **that operation had no
gene-content gate**. Every existing guard sat elsewhere: the 12F BUSCO trial gates
individual pool *replacements*, which on that run accounted for one upgrade and
ten additions, while the merge turned **972 backbone contigs into 468**. The 12H
purge gate covers purge_dups, which full mode never runs. The 12L do-no-harm
check compares size and telomere counts, both of which the merge improves.

- **New pass at 12G3**, immediately after the merge and before purge_dups. It
  compares complete BUSCOs in the backbone against the merge, identifies the
  backbone contigs carrying whatever was lost, and appends them with a
  `_gene_rescue` suffix. Validated against the real run: it identifies exactly
  the 112 lost genes and the 22 backbone contigs that carry them, 17,480,793 bp.
- **Bounded and honest.** A merge that loses nothing restores nothing. Restored
  contigs may duplicate sequence the pool already contributed — that is the
  intended trade, because a duplicated locus is recoverable downstream and a
  missing one is not. Every restoration is logged with the gene count, the point
  drop, and the base count. `STEP12_SKIP_MERGE_GENE_RECOVERY=1` disables it.
- **New `_busco_complete_by_id`** parses `full_table.tsv`, which lists one line
  per *copy*, so a duplicated gene correctly contributes every contig holding a
  copy.

A test now asserts that all three content-removing stages — the merge, purge_dups,
and the final comparison — each check gene content, so a future stage cannot be
added without one.

### Known limitation

The recovery pass restores whole contigs, so it can reintroduce redundancy: on the
run above it would add 17.5 Mb, taking the full deliverable from 256.5 Mb to about
274 Mb, or 1.08x the dikaryotic reference. Restoring only the gene-bearing
interval rather than the whole contig would be tighter and is not implemented.

---

## [1.5.0] — 2026-08-18

### `--assembly-mode {both,primary,full}` — both is the default

TACO now delivers **two genomes from one set of assemblies**: a nonredundant
primary representation in `mode_primary/` and a full sequence representation that
retains divergent alternate/haplotype sequence in `mode_full/`, plus
`final_results/assembly_modes_comparison.tsv` side by side.

`both` is the default because the two objectives disagree about what a good
assembly *is*, and nothing in the reads says which one a given project needs.
Choosing silently was the old behavior, and it was wrong in one direction for
every heterozygous, diploid, or dikaryotic sample.

Steps 0–10 describe the reads and the assemblers, so they are identical
whichever representation is delivered and run **once**. Only steps 11–14 depend
on the objective and run twice, each in its own working directory with the shared
step 0–10 products symlinked in — so the second deliverable costs no reassembly
and no duplicated 48 GB of FASTA. `--assembly-mode primary` or `full` emits just
one.

`full` is **not** a phased assembly. TACO does not order or orient contigs into
haplotypes and does not validate haplotype consistency. It retains alternate
sequence; it does not assign it.

`both` names no scoring profile and answers no removal question: the two
objectives disagree, so `AssemblyPolicy.sub_policies()` hands out one policy per
deliverable and asking a `both` policy whether purge_dups may run raises rather
than guessing.

### Full mode now has candidates that actually retain alternate sequence

`--assembly-mode full` selected among assemblies that had **all been reduced to
one representation before it ever saw them**. Several assemblers emit more than
one contig set and TACO used only the primary:

| assembler | also wrote | TACO used | now |
|---|---|---|---|
| hifiasm | `bp.hap1.p_ctg.gfa`, `bp.hap2.p_ctg.gfa` | `bp.p_ctg.gfa` only | `hifiasm_hap` = hap1 + hap2 |
| IPA | `a_ctg.fasta` (associate contigs) | `final.p_ctg.fasta` only | `ipa_alt` = primary + associate |

The consequence is visible in the real data: hifiasm and IPA reported 5.2% and
6.9% BUSCO duplication while canu and LJA reported 91.2% and 95.9%. The
full-representation winner was therefore decided by **which assemblers happen to
collapse haplotypes in their default configuration**, not by which recovers
alternate sequence best — and hifiasm and IPA had already produced exactly the
sequence the mode wants, then had it thrown away.

Both derived candidates cost **no extra assembly**: the files exist as soon as
the parent finishes. They are scored and selectable like any other candidate.

**They do not get a second vote.** The cross-assembler concordance check assumes
one ballot per independent assembly, and a derived candidate shares a graph, a
read set and an error profile with its parent — counting both would inflate
agreement and weaken the chimera evidence the vote exists to provide. The parent
votes; the derived candidate is scored but excluded from voting, and the run log
says which candidates were dropped and why.

Flye's `--keep-haplotypes` is a genuine third case but it changes the assembly
rather than reading a file that already exists, so it needs a reassembly and is
not enabled here. Assemblies produced before this release have no derived
candidates; re-run steps 4, 6 and 10 to generate them.

### purge_dups can no longer remove complete genes unnoticed

Running v1.5.0 on real data exposed a gap that predates it. purge_dups removes
haplotigs by sequence similarity and read depth and has no notion of gene
content, so it can remove the only copy of a locus — and **nothing downstream
would notice**. The 12F BUSCO trial gates *pool replacements*, not purging. The
12L do-no-harm check compares assembly size and telomere counts, and its size
rule explicitly *approves* a shrink that moves toward the declared genome size,
which is exactly what a purge does.

On *Puccinia triticina*, purging a peregrine backbone removed 39 haplotigs and
25.7 Mb, took BUSCO complete from **95.6% to 92.5%**, raised the T2T count from 4
to 6, and reported "final assembly quality OK".

This matters more from v1.5.0 on, and the release is what made it visible.
Backbone selection now tolerates a duplicated backbone *because purge_dups is
expected to clean it up*, which deliberately routes more sequence through the
unguarded step. The premise that duplication is "fixable, so don't pay for it at
selection" only holds if the fixing step is safe.

- **New gate at 12H.** purge_dups is accepted only if BUSCO complete falls by no
  more than the taxon tolerance already used by the rescue trial — 2.0 points for
  fungi, 4.0 for plants, 3.0 for vertebrates, 2.5 otherwise. On the run above
  that rejects the purge. A rejected purge keeps the unpurged assembly; the
  purged version and the removed haplotigs are preserved for review, and
  `STEP12_SKIP_PURGE_BUSCO_GATE=1` forces the old behavior.
- **The threshold table now lives in one place** (`_busco_delta_thresholds`), so
  a value tuned for the rescue trial cannot silently disagree with the purge gate.
- **12L reports gene content.** It reads the summary BUSCO already wrote, so the
  comparison costs nothing, and it warns when the drop exceeds the same tolerance.

The honest summary of the v1.5.0 primary genome on this data set is therefore
*better contiguity and consensus accuracy, worse gene completeness* — 35 contigs
against 117 and QV 76.8 against 47.9, but BUSCO complete 92.5% against 96.2%.
With this gate the purge that caused the loss is rejected, so the trade-off
should not recur; the run that produced those numbers predates the gate.

### Robustness fixes found by running it

**An external command could hang a run forever.** `run_cmd` had no time limit,
and the optional MUMmer `dnadiff` summary in the compare report is the one place
that matters: its `delta-filter -1` stage solves a 1-to-1 alignment selection
whose cost explodes with the number of alignments, and a repeat-rich genome
compared against a diverged assembly produces a great many. On a real run a
959 MB delta file left `delta-filter` burning a full core for 19 hours with a
zero-byte output, holding up step 14 behind an *optional* table. Two further
orphans from earlier runs on the same machine had been going for 8 and 11 days.

`run_cmd` now takes a `timeout`, and on expiry signals the command's whole
**process group** rather than just the immediate child — a shell pipeline's real
worker is a grandchild, which is exactly how those orphans outlived the runs that
started them. A timed-out command reports exit 124, as coreutils `timeout` does.
`dnadiff` is bounded at one hour, overridable with `TACO_DNADIFF_TIMEOUT`, and
the compare report proceeds without it.

**Deliverable-specific products were shared between representations.** The
per-mode workspace linked the contents of `busco/` and `merqury/`, which included
`busco/final`, `merqury/final` and `assemblies/final.telo.fasta` — the output of
whichever representation produced them. Step 13's BUSCO refresh hit
`Cannot call rmtree on a symbolic link`; worse, had it not failed, both
deliverables would have written their final QC through the same links and the
comparison table would have reported one representation's metrics as the other's.
Products under the `final` label are now excluded from sharing so each
representation builds its own, and both `rmtree` sites drop a symlink instead of
failing on it.

**A failure in one deliverable canceled the other.** `_run_step_list` exits
rather than raising on a failed step, and `SystemExit` does not derive from
`Exception`, so a crash in the primary representation aborted the whole run and
the full one never started — costing a day of compute on a real run. The
representations are independent products: each is now attempted, failures are
reported per representation, the surviving deliverable is still written to the
combined table, and the run exits non-zero.

**A timeout could crash instead of being handled.** The process-group id was
re-read after signaling, racing with the process exiting and raising an uncaught
`ProcessLookupError` on the very path meant to recover. It is read once.

### One combined result table at the run root

`final_results/final_result.csv` now carries the per-assembler columns from step
10 beside one `merged_<mode>` column per deliverable, so the nine assemblers, the
`--compare` reference, and both genomes read off a single table.
`assembly_modes_comparison.tsv` remains as the compact headline side-by-side.

This fixes a real defect in the first cut of `both`: every report was written
inside `mode_*/` and the root `final_results/` was left holding a **previous
run's** numbers, with nothing to say so. On a multi-day run that stale
`final_result.csv` reads as current for as long as the run takes. TACO now writes
`final_results/RESULTS_NOT_READY.txt` when the per-mode phase begins, naming the
superseded files and pointing at the `mode_*/` directories, and removes it once
the combined report exists.

### Backbone selection now scores what refinement cannot fix

The backbone is chosen *before* refinement runs. Measured on the real *Puccinia
triticina* run, backbone (ipa) → delivered:

| metric | backbone | delivered | |
|---|---|---|---|
| BUSCO C (%) | 96.3 | 96.2 | **unfixable** — only degraded |
| Merqury completeness | 91.74 | 91.51 | **unfixable** — only degraded |
| N50 (Mb) | 2.28 | 2.48 | barely movable |
| Merqury QV | 40.50 | 47.85 | fixable — polishing |
| BUSCO D (%) | 6.9 | 5.3 | fixable — purge_dups |
| strict T2T | 1 | 4 | fixable — telomere-pool rescue |

The pre-1.5.0 weights had this backward. `BUSCO_D` was charged 600/point — the
heaviest term after BUSCO itself — for duplication that step 12H removes anyway,
while contiguity got 30/contig and consensus accuracy 20 per QV point. On the
real data that rejected a **99-contig / N50 4.08 Mb / QV 71.7** candidate in
favor of a **329-contig / N50 2.28 Mb / QV 40.5** one, and delivered an assembly
*less contiguous and 24 QV points worse* than the candidate it passed over.

Selection now rewards BUSCO **complete** (not single — single is partly a
function of duplication, so rewarding it pre-purge rewards work already done),
weights contiguity and consensus accuracy far higher, and discounts duplication
and telomere counts because refinement supplies them. A 31-point QV spread is a
~1300× difference in consensus error, which 20/point effectively ignored.

**Merqury completeness means opposite things in the two modes**, which is why it
is not simply "a quality metric". It is the fraction of read k-mers present in
the assembly, and a collapsed haploid representation of a heterozygous sample
*cannot* contain both haplotypes'. So it is the dominant term in `full`, where it
measures the objective directly, and a light one in `primary`, where a value
above the haploid ceiling is evidence of retained redundancy rather than of
quality. On the real data the collapsed candidates sit at 91.7–93.8% and the
haplotype-retaining ones at 98.7–99.2%, so this distinction decides each mode's
choice. Without a dominant completeness term, `full` has **no** positive signal
for retention at all — BUSCO C saturates near 96% for every serious candidate and
duplication and size are deliberately unscored — and it selects a collapsed
147 Mb assembly over a genuine 254 Mb one.

Result on the real candidates: `primary` selects peregrine (99 contigs, N50
4.08 Mb, QV 71.7) and `full` selects lja (254.0 Mb, **1.002×** the published
253.5 Mb dikaryotic reference, duplication 95.9% against its 95.8%). The
46,812-contig artifact remains last in both.

**This deliberately changes primary-mode selection.** The pre-1.5.0 weights are
kept in `policy.LEGACY_PRIMARY_WEIGHTS` as a record of what the published v1.4.x
runs used, and a test asserts no profile silently resurrects them. Taxon presets
keep their relative ordering (fungal strictest on duplication, plant most
permissive because polyploidy inflates it) rescaled to the new level.

---

## [1.5.0-pre] — 2026-08-17

### `--assembly-mode {primary,full}`

TACO had one objective: a nonredundant primary representation. Every part of the
pipeline assumed it — the score rewarded `BUSCO_S` and penalized `BUSCO_D`, the
size term expected the haploid value, and purge_dups, self-dedup, and containment
filtering all removed sequence for looking redundant. On a heterozygous, diploid,
or dikaryotic sample that assumption is wrong in the same direction at every step,
and no single flag could undo it.

This release makes the objective explicit and selectable. **`primary` is the
default and reproduces pre-1.5.0 behavior exactly** — verified against the v1.4.1
formula across 72,000 randomized metric/taxon/genome-size combinations, maximum
difference 9.3 × 10⁻¹⁰ (float summation order).

`full` is **not** a phased assembly. TACO does not order or orient contigs into
haplotypes and does not validate haplotype consistency. The mode retains alternate
sequence; it does not assign it.

**Scoring** (`taco/policy.py`, new): the weights, the size expectation, and the
removal rules live in one module as pure functions of metrics and configuration,
so the two objectives cannot drift apart between `steps.py`, `backbone.py`, and
the reports. `backbone.py` carried a second, divergent copy of the formula; it now
delegates to the same function.

- `full` rewards `BUSCO_C` (S + D) instead of `BUSCO_S`, sets the duplication
  penalty to **zero** — neutral, not rewarded — and doubles the Merqury
  completeness weight to 400. QV stays separate from completeness at 20 in both
  modes.
- Taxon overrides can no longer defeat the mode: `--taxon fungal` would otherwise
  reimpose `w_busco_d = 600` on `full`.
- `BUSCO_C` is now read from `assembly_info.csv`, and reconstructed as S + D on
  older tables that lack it.

**Size no longer ranks in `full` mode.** No assembler declares an expected total
size for a heterozygous sample — hifiasm and Verkko take no genome size, Flye's is
optional, and Canu's `genomeSize` is a coverage knob feeding `corOutCoverage` and
`mhapSensitivity`, not a size expectation. Where a declared size is compared
against a finished assembly, NCBI uses a wide taxon-derived band anchored on the
haploid value to flag for human review, never to rank. `full` mode reports the
band `[-g, 2 × -g × 1.15]` and scores nothing from it.

The over-assembly guard is structural instead. On a real *Puccinia triticina* read
set (`-g 127m`, `--taxon fungal`), a 287 Mb / 46,812-contig artifact had the
highest Merqury completeness of any candidate (99.26 %) and sat *inside* any band
wide enough to admit the genuine 254 Mb assembly, so size could not separate them.
Contig **density** does: 163 contigs/Mb against 3.8 and 2.5. Density rather than
absolute count, because an assembly carrying twice the sequence carries twice the
contigs, and charging the absolute count penalizes it for exactly the property the
mode selects for. Fragmentation is priced at the same per-contig rate in both
modes — 30 per contig, evaluated against a haploid-sized assembly — so only the
double-charge is removed, not the guard.

**Selection alone is not enough.** `full` mode also declines every step that
removes sequence *because it looks redundant*: purge_dups (12H), post-upgrade
dedup (12G2), Tier 2 self-dedup, and containment filtering of the telomere pools
(step 11). Each logs why it was skipped.

Contaminant screening (13A) and chimera resolution (13B) run in **both** modes.
Foreign sequence and read-unsupported joins are errors under either objective and
are identified from composition, read depth, and read spanning — not from a contig
resembling a second haplotype.

`--no-purge-dups` is unchanged and still supported. Precedence: it can only ever
*disable* purge_dups, never enable it. `full` mode already skips the step, and
omitting the flag does not bring it back.

**Reports.** New `final_results/assembly_mode_summary.txt` states
`Assembly representation mode: primary|full` plus the selected assembler,
selection score, active score profile, BUSCO C/S/D, Merqury QV and completeness,
assembly length, contigs, N50, T2T and single-end telomere counts, and the
individual score contributions behind the winning score. The mode also appears in
`final_result.csv` and `selection_decision.txt`, and
`assemblies/selection_debug.tsv` gains a `score_*` column per scoring component,
so a run records *why* an assembly won rather than only that it did.

**Advisory.** A `primary` run that ends with Merqury completeness ≥ 95 %, BUSCO
duplication ≥ 20 %, and assembly size ≥ 1.4 × the expected haploid value emits a
warning suggesting `--assembly-mode full`. It is worded as a prompt to look, not a
biological conclusion: these metrics alone do not establish the ploidy of the
sample, and the warning does not claim they do.

`smart` and `both` modes are designed for but deliberately not exposed yet; a test
asserts the CLI rejects them.

41 new tests in `tests/test_assembly_mode.py` (133 total, up from 92), covering
CLI defaults and rejection, v1.4.1 score equivalence, the two modes disagreeing on
one input, the real *P. triticina* selections in both modes, the removal policy,
and the advisory's wording.

### Scientific-validity fixes found by re-analyzing the real data

Three problems surfaced while validating the release against the *Puccinia
triticina* Pt76 run, and are fixed here.

**The telomere reward was unbounded in a quantity biology bounds.** A contig has
two ends, so it carries at most two telomeres, and a genome holds
2 x n_chromosomes x ploidy of them. The score rewarded `T2T` and single-end
telomere counts linearly with no ceiling. On the real data the MBG assembly
reported **1,668** single-end-strong telomere contigs for an 18-chromosome genome
where at most 36 are possible (72 as a dikaryon), earning +250,200 points —
three times its entire BUSCO term. Only its 46,812 contigs kept it from winning,
so the safety margin rested entirely on fragmentation rather than on any bound.

`policy.telomere_reward_cap` now bounds the rewarded count at
2 x (expected total / `MIN_PLAUSIBLE_CHROMOSOME_BP`), with the floor set to 1 Mb.
That is deliberately permissive — it implies 127 chromosomes in a 127 Mb haploid
genome, roughly seven times the real 18 — so the bound only fires on counts that
are impossible by a wide margin. Verified to change **no** real selection: both
modes rank the nine real candidates exactly as before, and only the artifact's
score moves. Counts past the bound stop being rewarded rather than being
penalized, and with no `-g` no bound can be derived so behavior is unchanged.

This is the one place primary mode is no longer identical to v1.4.1. The
equivalence test now asserts equality within the bound and tests the bound
separately, because past it v1.4.1 was rewarding false positives. A tighter
taxon-aware bound would cut deeper but would start capping legitimate
candidates, so it is left as future work.

**Contaminant screening could misfire on a haplotype-retaining assembly.** Such
an assembly is inherently bimodal in depth: homozygous regions collapse to one
contig at full read depth, heterozygous regions are kept twice at about half. If
the collapsed population is the heavier cluster it becomes the core, and a
legitimate half-depth haplotype contig then sits at ratio ~0.5 against a 0.34
removal threshold — a 1.5x margin. Full mode now screens with
`DEPTH_RATIO_MAX_FULL = 0.25`, restoring a ~2x margin. Screening is still never
skipped, and the change can only make removal more conservative.

**The compare report read representation differences as mis-assemblies.** When a
collapsed assembly is compared against a reference that keeps both haplotypes,
each reference haplotype pair maps onto the same single contig, and
`contig_to_contig.tsv` labeled those rows `1-to-N (split)` — which reads as
though TACO broke a chromosome. Against the Pt76 reference (253.5 Mb, 36
chromosomes = 2 x 18), **97%** of chromosome-scale reference contigs shared a
best target and 48 rows were labeled "split". New
`final_results/compare_report/REPRESENTATION_NOTES.txt` reports the shared-target
fraction and, above 50%, explains that these are two representations of one
genome, warns that `weak_regions_compare.tsv` and `unique_compare_contigs.tsv`
will list sequence absent by design, and says so in the run log. It makes no
defect claim about either assembly.

### Known issue

`taco/steps.py:_filter_redundant_to_protected` is defined but never called. It
predates this release and was left in place.

---

## [1.4.1] — 2026-08-03

### Audit only: the Tier 1 / cross-assembler disagreement is now recorded

Step 12 assigns Tier 1 immutability from telomere status alone and never
consults the cross-assembler verdict. A contig fusing two chromosomes at ends
that are themselves incomplete carries genuine telomeres at both outer ends, so
it classifies strict-T2T, becomes immutable, and the telomere pool is barred
from competing for that chromosome — precisely because the contig fakes the
property that earns protection. 13B resolves it later, but refinement has
finished by then and no replacement can be offered for the pieces.

This release makes that visible and changes nothing else.

- New `assemblies/tier_assignment.tsv` (copied to `final_results/` when present)
  with `tier` — what TACO used — beside `would_be_tier` — what the verdict
  implies — plus the vote counts and the reason. Where the two columns differ,
  the run directory now says so.
- A warning naming each contested Tier 1 contig and stating that resolution is
  deferred to 13B.
- Reuses `_purify_all_contig_concordance`, which already excludes the backbone's
  own assembly from the voter set and already handles a missing minimap2 or
  fewer than two voters. Cost is one call, roughly 25 s on a 45 Mb fungal
  genome; skipped entirely under `--concordance-mode off`.

### Correctness fixes from a v1.4.0 audit

An audit of v1.4.0 produced 70 candidate findings, of which 15 were confirmed by
adversarial re-reading of the source. Six are fixed here; the rest are recorded
in the audit report and deferred.

- **Consensus breakpoint was computed from every voter, not only the ones that
  saw a split.** Assemblies voting `intact` or `uninformative` contribute no
  interior gap but still entered the estimate, dragging the consensus toward the
  middle of the contig and cutting at a coordinate no assembly proposed. The
  estimate now uses only `split` voters. This is the most consequential fix in
  the release: it changes where a chimera is cut.
- **A contig could be removed on read depth alone.** Removal requires two
  independent primary signals, or one plus corroboration — but the only live
  corroborator, the fraction of the contig at near-zero depth, is the left tail
  of the same depth distribution whose median is the primary. Both are depressed
  together whenever reads simply fail to map, so a divergent haplotype at low
  coverage with entirely normal composition was classified `foreign` and dropped
  from the delivered assembly. Coverage uniformity may no longer corroborate a
  depth primary. The regression test that should have caught this was passing
  against a fixture without `zero_bp`/`low_bp`, a shape the pipeline never
  produces; it now carries them.
- **The MAD was count-weighted around a length-weighted center.** A numerous
  small-contig population could inflate the depth sigma past `centre / 3.5`,
  making the depth signal unable to fire for any contig in the assembly. The
  scale is now weighted exactly as the center is.
- **`purify_excluded.fasta` was written only when something was removed** and
  never cleared, so a re-run that removed nothing left the previous run's file
  claiming exclusions from an assembly that still contained them. It is now
  rewritten unconditionally, empty when nothing was removed.
- **A single donor could be consumed by several backbone targets**, deleting N
  contigs and inserting the same donor sequence N times. A donor is now used at
  most once.
- **`chimera_decision` labeled every non-mis-join verdict `corroborated`**,
  asserting cross-assembler agreement for contigs that were never assessed.
  Only a genuine corroborated verdict now says so.

Tests: 92 across four files (five added). Each new test was verified to fail
against the unfixed code.

**No output changes.** No tier is modified, no path is altered, and the
delivered assembly is byte-identical to v1.4.0 on the same inputs.

Acting on the disagreement is deferred deliberately. Demoting a contested contig
to Tier 2 makes it eligible for donor paths — 12E2 single-telomere rescue at
`RESCUE_MIN_COV_BB` 0.60, and the 12F2 partial path at `min_tcov` 0.15 — whose
thresholds have not been evaluated against a multi-chromosome target. Demotion
is therefore not obviously safer than current behavior, and needs a genome where
a donor can actually win before it is worth the risk.

---

## [1.4.0] — 2026-07-31

A minor-version bump rather than a patch, because the default behavior changes in
two ways a user will notice: `final_assembly.fasta` is now the purified assembly
rather than the unfiltered merge, and step 13 modifies the assembly where it
previously only measured it. Both are opt-out via `--purify-mode off`. The work
was developed on `main` as 1.3.9 and released as 1.4.0; there is no 1.3.8 or
1.3.9 release.

### Purification now runs before measurement, and it acts

v1.3.7 detected both of the problems below accurately and then changed nothing a
user reads. The contaminant screen ran in step 14A, *after* step 13 had already
measured BUSCO, telomere counts, QUAST and Merqury on the unfiltered merge, and
`final_assembly.fasta` was copied from that same unfiltered merge. On the run
that motivated it the delivered assembly was 50,873,779 bp at 49.34% GC when the
truth was 46,238,257 bp at 47.62% — the screen had identified the 4.64 Mb
responsible and written it to a side file nothing downstream consumed.

- **Step 13 is now "Purify + final QC"**, with sub-steps 13A screen, 13B chimera
  resolution, 13C emit, 13D measure. Purification and measurement are the same
  step, so the reported metrics cannot describe a different assembly than the one
  TACO delivers — the v1.3.7 defect is now structurally impossible rather than
  merely fixed. **Step numbers are unchanged**: every existing `-s` recipe keeps
  working and silently becomes correct.
- **`final_assembly.fasta` is the purified assembly.** `final.merged.fasta`
  remains the unfiltered merge and everything removed is preserved in
  `purify_excluded.fasta`, so no sequence is destroyed.
- New `taco/purify.py` holds the decision logic as pure functions, unit-tested
  without minimap2, samtools, BUSCO, or a reference genome.

### Contaminant screening rewritten to generalize beyond one genome

The v1.3.7 rule — fixed `DEPTH_RATIO_MAX = 0.34` **or** `GC_DEV_MAX = 8.0` on a
length-weighted median baseline — worked on the genome it was written for and is
unsafe elsewhere. Three defects, each measured:

- **The baseline inverted.** A length-weighted median is the assembly's midpoint
  by length, so once foreign sequence passes half the assembly the "baseline"
  becomes the contaminant and the host is measured against it. Iterative trimming
  cannot rescue this — the initial center sits on the wrong cluster and the sigma
  is inflated by the very gap that should be resolved, so nothing is ever trimmed.
  The baseline is now the **heaviest cluster** by sequence length, which is also
  what "modal depth" means in a coverage histogram.
- **A single signal could condemn a contig.** No established tool removes sequence
  on GC or coverage alone; BlobToolKit's premise is that GC × coverage needs a
  third axis because both are confounded. Removal now requires **two independent
  primary signals**, or one plus corroboration.
- **Fixed absolute thresholds do not travel.** 8 GC points is ~13 robust sigma on
  this genome and under 1 sigma on a compositionally heterogeneous one. Rules are
  now a robust z-score **and** an absolute effect size — both, because either
  alone fails. A textbook Iglewicz–Hoaglin `|z| > 3.5` on GC flags a real 3.08 Mb
  chromosome here at −3.71; the 5-point floor is what saves it. Robust sigmas are
  floored, because the real contigs of a good assembly agree on depth to within
  half a percent and the unfloored MAD collapses.

Guards that protect real-but-atypical sequence, each of which was a live false
positive on the motivating genome: telomere arrays at both ends veto a removal;
depth ≥ 3× core marks an organelle or collapsed repeat rather than a contaminant;
composition-only calls need ≥ 250 kb; and nothing is removed if it would discard
more than 25% of the assembly or if no cluster holds a majority of the length.

**BUSCO gene content is recorded and deliberately never used to remove anything.**
On this genome a real 3.08 Mb chromosome carries zero ascomycota BUSCOs while the
bacterial contig carries three, so gene content does not separate them. The
v1.3.7 `busco_count == 0 and length >= 1 Mb` signal would have condemned a
chromosome.

- New `--metagenome`: report anomalies, never remove. A metagenome has no single
  modal depth or composition for the screen to work against, so the premise does
  not hold. For co-cultures, holobionts, lichens, and symbiont assemblies.
- New `--purify-mode {on,off}`, default `on`. `--no-contam-screen` is kept as an
  alias for `off` so a v1.3.7 command line does not silently gain a stage that
  modifies the assembly.

### Chimeras are now detected across every contig, and repaired

v1.3.7 cross-checked only contigs already classified strict-T2T, purely to correct
the backbone selection score, and deleted the alignment PAFs immediately after
voting. So a chimera without telomeres at both ends was never examined, the
coordinates needed to act were discarded, and the corrected count reached the
*score* but never the *reported* metric — `final.telo_metrics.tsv` still claimed
9 strict-T2T contigs when 8 was true.

Step 12C's two chimera gates could not have caught it either, for structural
reasons rather than threshold ones:

- It screened only `protected_telomere_contigs.fasta`; the offending contig
  arrived by the **backbone** path, which 12C never examines.
- Its cross-assembly mapping gate treated the backbone's *own* assembly as a
  voting target, so every backbone contig self-aligned at 100% coverage. The gate
  could not fire on any assembler-sourced contig.
- Its size gate derived the threshold from the largest contig across assemblies —
  itself another assembler's mis-join here, inflating the threshold to 13.0 Mb.
  The chimera was in any case *smaller* than the largest real chromosome, so no
  length statistic can catch it.

13B replaces all of that:

- **Every** contig above 500 kb of the delivered assembly is cross-checked
  against every other assembly, excluding the backbone's own.
- Per-voter alignment intervals are **retained**, which is what makes a
  reference-free breakpoint estimate possible: each voter that splits the contig
  into two substantial, largely disjoint components contributes the midpoint of
  its interior gap, and the consensus is the median of those.
- **Read-spanning evidence decides.** Coverage continuity is not evidence of a
  real join — Tigmint documents exactly this case — and at a repeat-mediated
  mis-join the depth *spikes* rather than dips, so a coverage-dip test sees
  nothing, while reads map inside the repeat rather than clipping at its edge, so
  a clip-pileup test sees nothing either. The spanning count is judged as a
  fraction of the local pileup, because depth varies by an order of magnitude
  across such regions.
- Nothing is broken without a cross-assembler majority, an agreed breakpoint,
  **and** reads that fail to span it. If reads *do* span it, TACO keeps a join
  every other assembler missed. Telomeres at both ends raise the bar to a read
  refutation rather than exempting the contig, since a fused pair of incomplete
  chromosome ends is precisely what produces that signature.
- New `--chimera-action {split,replace,report,off}`, default `split` — cutting in
  place preserves the backbone sequence and its polish. New `--spanning-anchor`.

### Also fixed

- `concordance_verdict` treated a tie as corroboration (`n_split > n_intact`), so
  an even split silently protected a candidate. A tie now flags.
- `split_vote` called "split" on two targets piled onto the *same* part of the
  query, contradicting its own docstring. Overlapping components now read as a
  repeat family, not a junction.
- GC is computed over ACGT only. Counting N as non-GC dragged a gappy contig's GC
  downward and could manufacture a composition outlier.
- `taco/concordance.py` now keeps only the voting logic.
- The 1.3.7 entry below undercounted its own tests (26 in `test_concordance.py`,
  37 across the suite, and it omitted `test_v135_rescue.py`).

### Verified against the reference rather than asserted

Sub-steps 13A–13C were run against the real *Fusarium tricinctum* MsR-QD66
assembly (HiFi SRR33612568) and checked against the published Hi-C assembly
GCA_050859235.1, which purification never sees:

| | v1.3.7 delivered | v1.3.9 delivered | Hi-C reference |
|---|---|---|---|
| contigs | 11 | 10 | 10 |
| total length | 50,873,779 | 46,238,257 | 46,236,445 |
| GC | 49.34% | 47.62% | 47.61% |
| strict T2T reported | 9 | 8 | 6 |

The two removed contigs are exactly the two QUAST reports as wholly unaligned to
the reference (4,635,522 bp, matching to the base pair). The chimera was cut at
1,593,380, inside the reference-derived seam at 1,577,246–1,598,296, on 7-of-7
assembler agreement with 0 of 6,040 overlapping reads spanning it — while the same
test *supported* the equivalent position on all eight other contigs, at 34 to 295
spanning reads each. The honest telomere count of 8 still beats the published
Hi-C assembly's 6.

The mechanism, for the record: the seam is a collapsed tandem array of a ~7.9 kb
unit at ~19× the assembly's modal depth, spanning roughly 62 kb as assembled and
therefore ~1.2 Mb in reality. Peregrine joined two chromosomes through it. Both
reference scaffolds terminate at that array, which is why each lacks exactly the
telomere facing the other.

### Provenance

A **Use Of AI Assistance** section is added to the README, noting that an AI
coding assistant was used for bug fixing and for implementing the chimera check
and contaminant screening, and that no sequencing data, assembly or reference
genome was generated or modified by an AI tool.

Tests: `tests/test_purify.py` adds 58, including the three false-positive traps
above as regressions, and `tests/run_all.py` runs the whole suite in one command for the first time. Suite total 87 across four files (eight contaminant tests moved out of `test_concordance.py` into the rewritten module's own file), no external tools required.

---

## [1.3.7] — 2026-07-27

### Cross-assembler concordance: a single assembler's mis-join can no longer win backbone selection

- **New `taco/concordance.py`.** A contig that one assembler presents as
  telomere-to-telomere while every other assembler splits it in two is a chimera
  candidate. This matters for backbone selection because fusing two *incomplete*
  chromosome ends yields a contig with genuine telomeres at both outer ends and
  nothing anomalous at the seam, so it satisfies the strict-T2T test and earns the
  T2T score reward. The join is only visible as disagreement between assemblers —
  evidence a single-assembly tool does not have.
- **`_t2t_concordance_check` (step 10).** Aligns every strict-T2T contig against
  the other assemblies and takes a majority vote, writing
  `assemblies/t2t_concordance.tsv` and `assemblies/t2t_corroborated.tsv`.
  `_auto_select_backbone` now scores on the *corroborated* T2T count. Requires no
  reference genome. New `--concordance-mode {exclude,flag,off}` (default
  `exclude`). Two split votes are required and uninformative voters are ignored,
  so a single fragmented assembly cannot flag a real chromosome; if alignment
  fails the verdict is `unresolved` and no contig is penalized.
- **Mis-join detection in compare mode.** `contig_to_contig.tsv` recorded a contig
  absorbing two compare chromosomes as two independent "1-to-1" rows, making the
  fusion invisible. `mis_join_candidates.tsv` now reports it with the implied
  junction coordinate.
- **Contaminant screen.** Flags final-assembly contigs whose read depth and/or GC
  mark them as foreign and writes `contamination_candidates.tsv` plus
  `final.clean.fasta` **alongside** the unfiltered assembly. purge_dups removes
  haplotypic duplication and will not remove a foreign genome. New
  `--no-contam-screen`.
- **Tests restored and extended.** `tests/` was deleted in 078f4ae; restored
  `test_review_fixes.py` (8 tests) and added `test_concordance.py` (24 tests).
  32 tests, no external tools required.

Found in a real run: on *Fusarium tricinctum* MsR-QD66 (SRR33612568), Peregrine
joined reference chromosomes 6 and 10 into one 6.5 Mb "T2T" contig that seven
other assemblers and the published Hi-C assembly (GCA_050859235.1) keep separate.
The concordance check flags it (4/4 voters split) and corrects the count 9 -> 8
without changing which assembly is selected. The contaminant screen flags 4.64 Mb
of probable bacterial sequence (66.6% GC, 19x against 312x modal depth, three 16S
rRNA copies) that purge_dups had retained.

---

## [1.3.6] — 2026-07-20 — correctness & robustness review

Fixes from a full-pipeline correctness / scientific-validity / robustness
audit. No CLI or output-schema changes; behavior is more correct.

### Scientific validity

- **Telomere end scoring is now 5'/3' symmetric.** `score_end` measured the
  distance-to-terminus from the window's left edge for *both* ends, so a
  telomere array shorter than the score window scored lower at the 3' end than
  the identical array at the 5' end. This systematically under-called
  right-terminal (and therefore T2T) contigs — worst for fungi (short arrays,
  300 bp window). The right end now measures distance from the right edge.
- **BUSCO metrics are gene-based again.** The `full_table.tsv` parser counted
  one line per *copy*, inflating Duplicated and the denominator `n` so all
  percentages diverged from BUSCO's own `short_summary`. Counts are now
  de-duplicated per BUSCO id. This feeds backbone selection (heavily weighted),
  so it could previously flip the chosen assembler.
- **Quickmerge chimera gate uses the two actual parent contigs**, not the
  largest contig in the whole file, so an inter-chromosomal chimera of two
  mid-sized contigs can no longer slip into the protected T2T pool. A validated
  merge now also requires *both* parents to contribute substantial coverage
  (was: either one), matching the documented "genuine two-parent join" rule.
- **Dedup no longer deletes both copies of an equal-length duplicate.**
  `_fasta_clean_contained` and `_self_dedup_non_telomeric` broke length ties
  such that both members of an equal-length pair were dropped; they now keep
  one representative and `_fasta_clean_contained` gained an identity gate and
  cross-block coverage aggregation so a shared repeat can't delete a distinct
  contig.

### Correctness / robustness

- **Protected Tier 1 contigs can no longer be silently overwritten** during
  step-12 rescue/upgrade: a donor whose name collides with a *different*
  backbone contig is inserted under a unique key instead of clobbering it.
- **Gzipped `--reference` / `--compare` downloads are decompressed** by gzip
  magic-byte detection (a `.fa.gz` URL was saved as `.fasta` and read as raw
  gzip). Download commands now use argv (no shell) so tokened URLs are safe.
- **purge_dups read alignment detects minimap2 failure** instead of masking it
  behind a `gzip` pipe (a truncated PAF previously biased the coverage cutoffs).
- **Merqury reuses only a database built at the resolved k**, never a stale
  db at a different k-mer size.
- **Polishers stream to disk** instead of buffering the whole genome in memory.
- **Compare-report variant calling** passes `--cs` so `paftools.js call`
  produces a real VCF instead of an always-empty one.
- **CLI:** an explicit `--telo-score-window 500` is honored (was overridden by
  the taxon default); `--assembly-only` combined with `--steps` now intersects
  and warns instead of silently discarding the selection; the standalone
  `telomere_detect` / `telomere_pool` CLIs apply their `--taxon` and
  threshold flags.
- **Docs/env:** `seqkit` (used by Step 0) added to `taco-env.yml`; unused
  `seqtk`/`bwa` removed from version logging and install docs; README version
  strings and the fungal BUSCO default (`fungi_odb10`) corrected; `setup.py`
  reads README as UTF-8; `utils.revcomp` complements the full IUPAC alphabet.
- Version bumped to 1.3.6 across `__init__`, `cli`, `pipeline`, `steps`,
  `setup.py`, and the README. New regression tests in
  `tests/test_review_fixes.py` lock in the fixes above.

---

## [1.3.5] — 2026-07-08

### Telomere-aware rescue: short terminal telomere caps are no longer discarded

- **`step_12` single-end rescue rejected donors that only add a short telomere
  cap.** `_screen_rescue_candidates` dropped any donor whose extension beyond
  the backbone was `< RESCUE_MIN_EXT` (1000 bp; `ext = max(qs, qlen-qe)`),
  regardless of whether that overhang was a telomere. Telomere repeat arrays are
  frequently only a few hundred bp (especially fungal), so a donor whose sole
  new contribution was the terminal telomere the backbone lacked was thrown out
  for being "too short." Symptom: single-end backbone contigs stayed uncapped
  even though a telomere-bearing contig from another assembler aligned cleanly to
  the terminus and extended it by ~100–575 bp (one run rejected 15 such donors
  solely for `ext<1000`).
  - `_parse_paf_rescue_hits` now computes, per hit, the strand-aware donor
    overhang at the touched backbone terminus and whether it carries a telomere
    (`term_ext`, `overhang_has_telomere`, `overhang_telo_score`), reusing TACO's
    own `score_end` scorer with the taxon's motif families (new
    `_overhang_is_telomeric` helper).
  - `_screen_rescue_candidates` exempts a short extension from the
    minimum-extension rejection when the overhang adds a telomere the backbone
    lacks, and ranks telomere-adding donors highest. All chimera guards
    (identity, aligned_bp, backbone/donor coverage, terminal touch) and the
    downstream BUSCO trial (12F) are unchanged, so accepted caps still cannot
    reduce completeness.
  - `single_tel.replaced.debug.tsv` and `single_tel.candidates.tsv` gained the
    telomere-overhang columns, so each rejected/accepted end is auditable.

### purge_dups over-purge guard (do-no-harm on genome size)

- **`step_12` purge_dups could shrink a well-sized assembly and still pass QC.**
  `_purge_dups_safety_check` was size-only, rejecting only on an absolute floor
  (`min_expected_ratio × genome size`) and a `max_bp_drop` fraction. A
  low-duplication backbone (BUSCO D ≈ 0) could therefore lose ~9 Mb / dozens of
  single-copy BUSCOs and be accepted (drop within the 25% limit, above the 85%
  floor). Added an **overshoot guard**: if purge pushes the assembly *further
  below* the expected genome size than it already deviated while dropping > 5%,
  the purge is rejected and the unpurged assembly is kept
  (`reason = overshoot_below_expected_size`). Legitimate purges that move an
  inflated assembly *toward* the expected size are unaffected. Override with
  `PURGE_DUPS_ALLOW_OVERSHOOT=1`.

### Notes

- Version bumped to 1.3.5 across `__init__`, `cli`, `pipeline`, `setup.py`, and
  step logs/reports. No change to public step numbering (0–14) or CLI flags.

---

## [1.3.4] — 2026-05-06

### Two major QC-reuse bugs (`--busco` override + Merqury cross-contamination)

- **`--busco <lineage>` was silently ignored on re-runs.**
  ``step_08_busco`` checked "is there ANY BUSCO output for this
  assembler?" but not whether the existing output was produced by the
  lineage the current run asks for.  Concrete symptom: a first run with
  ``--taxon fungal`` writes ``busco/canu/run_fungi_odb10/``; a second
  run with ``--busco ascomycota_odb10`` finds the stale dir, logs
  "Found existing BUSCO metrics for canu → using busco/canu", and
  reports fungi_odb10 numbers under the ascomycota_odb10 lineage row.
  Fixed by lineage-aware reuse detection: TACO now requires the
  on-disk output to come from the requested lineage (matched against
  ``run_<lineage>/`` and ``short_summary.*.<lineage>.<base>.{txt,json}``
  with content-based fallback for older BUSCO versions).  When the
  lineage doesn't match, the existing dir is wiped and BUSCO re-runs.
- **Merqury rows were all reusing canu's QV/completeness.**  In
  ``_merqury_metric_candidates`` an over-broad ``merqury/**/*.qv``
  recursive glob fired whenever the prefix's ``out_dir`` resolved to
  the merqury root (the legacy flat-prefix case), so every label
  beyond canu silently matched ``merqury/canu/canu.qv``.  Symptom in
  ``final_result.csv``: every assembler row (compare, flye,
  nextDenovo, raven, ...) reported the same 50.2131 / 98.6575 from
  canu's files.  Fixed in two ways: (a) the recursive ``**`` glob is
  now scoped to the per-assembler subdir (``out_prefix/`` only, never
  the merqury root); (b) a belt-and-suspenders filter rejects any
  candidate whose basename doesn't contain the label token, so
  ``canu.qv`` cannot match label "compare" / "flye" / "raven".

### Step 14C compare report — telomere-aware chimera diagnosis

- **Telomere status now flows into every layer of the compare report.**
  Step 14C reads both `assemblies/final.telomere_end_scores.tsv` and
  `assemblies/compare.telomere_end_scores.tsv` and propagates per-end
  telomere flags (`yes` / `no` / `n/a` for "not tested at this
  boundary") into:
  * `contig_to_contig.tsv` gains `compare_left_telo`,
    `compare_right_telo`, `compare_telo_tier`, `best_target_left_telo`,
    `best_target_right_telo`.
  * `contig_to_contig_pairs.tsv` gains four boundary-touch flags
    (`touches_compare_left/right`, `touches_final_left/right`) plus
    four pair-localized telomere flags
    (`compare_left_telo_at_pair`, `compare_right_telo_at_pair`,
    `final_left_telo_at_pair`, `final_right_telo_at_pair`). A flag of
    `n/a` means "the alignment block doesn't reach this boundary so
    the telomere status of that boundary is not relevant for *this*
    pair" — distinct from "no telomere".
  * `split_mappings.tsv` gains a `chimera_evidence` text column that
    summarizes the telomere pattern at the four boundary points and
    classifies the split as one of: full-coverage split, strong
    chimera signal (full chromosomes joined via shared telomeres),
    chimeric extension into a fragment, or split-likely-real
    (compare is incomplete on the side without a telomere).
- **Verified on the user's KMAF11 run** — `JALGOQ010000001.1` (6.51 Mb)
  is now correctly diagnosed as a chimeric extension: compare has
  telomeres at both outer ends, contig_2 (TACO's strict-T2T half)
  contributes a telomere at its inner end, contig_1 (TACO's
  single_tel_strong half) does not — exactly the "chrA fully joined
  to chrB-fragment" signature.

### Step 14C compare report — split-mapping detection (chimera-aware)

- **Bug fix: `contig_to_contig.tsv` previously hid 1-to-many mappings.**
  When a single `--compare` contig spanned more than one final contig
  (i.e. the published assembly chimerically joined two real chromosomes,
  or TACO split one chromosome where the published assembly didn't), the
  per-query reducer kept only the highest-coverage target. The row read
  as a clean `1-to-1` even though `aligned_bases` exceeded `target_len` —
  the tell-tale arithmetic of a hidden split. Reproduced on a real
  Venturiaceae run: `JALGOQ010000001.1` (6.51 Mb) silently mapped to
  `contig_1` (3.50 Mb), hiding its 3.05 Mb second half on `contig_2`.
- **Per-pair aggregation surfaces every significant compare→final
  partner.** Each row in the new `contig_to_contig_pairs.tsv` is one
  (compare_contig, final_contig) pair above the significance bar
  (max of 20 % of the compare contig OR 100 kb absolute), with
  per-pair compare coverage, final coverage, identity, and dominant
  strand. `JALGOQ010000001.1` now writes two rows: `contig_1` (53.1 %,
  strand `-`) and `contig_2` (46.9 %, strand `+`) — opposite strands
  is a strong inversion/chimera signal.
- **`contig_to_contig.tsv` gains three columns**: `n_significant_targets`,
  `relationship` (`1-to-1` / `1-to-N (split)`), and `secondary_targets`
  (semicolon-separated `name:bp:cov%` list of partners besides the
  best one).
- **`split_mappings.tsv` new file**: only compare contigs that map to
  more than one significant final contig. Inspect this first to triage
  candidate chimeric joins. On the Venturiaceae paper assembly, two
  rows appear: the 6.5 Mb chimera and a smaller 1.9 Mb split into
  contig_18 + contig_20.
- **`synteny_blocks.tsv` now reports four relationship classes**:
  `1-to-1`, `1-to-many (compare splits across N final contigs)`,
  `many-to-1 (M compare contigs → 1 final)`, and `many-to-many` —
  derived from the per-pair table so splits in either direction are
  visible.

### Step 14C compare report — unique-contig listings, synteny blocks, Circos config

- **`final_results/compare_report/unique_compare_contigs.tsv`** lists every
  contig in the `--compare` FASTA whose best alignment to `final.merged.fasta`
  covers less than 5 % of its length (or has no alignment at all).
  Companion file **`unique_final_contigs.tsv`** does the same for
  `final.merged.fasta` contigs that received < 5 % coverage from any
  compare alignment. Each row carries length, coverage, best partner,
  identity, and a `reason` column (`no_alignment` /
  `coverage_below_5pct`).
- **`synteny_blocks.tsv`** aggregates 1-to-1 best mappings (compare
  coverage ≥ 50 %) and tags each block as `1-to-1` or `many-to-1`,
  giving an at-a-glance synteny summary instead of forcing the user to
  derive it from the per-alignment PAF.
- **`compare_report/circos/`** is a self-contained, ready-to-render
  Circos input bundle. TACO writes:
  * `karyotype.txt` — every compare and final contig as an ideogram
    (`C_<name>` for compare contigs, `F_<name>` for final contigs).
  * `links.txt` — one ribbon per minimap2 alignment ≥ 5 kb.
  * `circos.conf` — minimal Circos configuration referencing the two
    files above; uses the standard Spectral palette so compare and final
    contigs are colored differently.
  * `README.txt` — instructions to install Circos and run
    `cd circos && circos -conf circos.conf`.
- **Optional `compare_paftools_variants.vcf`** — when `paftools.js` and
  `k8` are on PATH (they ship with minimap2), TACO sorts the PAF by
  target and runs `paftools.js call -L 10000` to emit assembly-vs-assembly
  SNVs and SVs in VCF format. Skipped silently when the binaries are
  missing.
- **Optional `compare_mash_distance.tsv`** — when `mash` is on PATH,
  TACO writes the scalar Mash distance between the two assemblies as a
  quick "are these the same organism" check.

### `--compare` and `--final-fa`: passive comparison + arbitrary final assembly

- **New `--compare <fasta>` flag.** A user-supplied FASTA can now be brought
  into the comparison tables (BUSCO, Telomere, QUAST, Merqury) without any
  risk of contaminating the refinement pipeline. `--compare` is excluded
  from backbone selection, the telomere-pool quickmerge candidate set,
  polishing, and `purge_dups`. Step 14 adds a new sub-step **14C** that
  runs whenever `--compare` is set: minimap2 `-cx asm5` aligns the compare
  FASTA to `final.merged.fasta` and TACO writes
  `final_results/compare_report/contig_to_contig.tsv` (per compare
  contig: best 1-to-1 target, aligned bases, identity %, query coverage)
  and `final_results/compare_report/weak_regions.tsv` (10 kb windows of
  `final.merged.fasta` whose compare coverage is below 50 %). When QUAST
  is available, `compare_quast/` adds reference-based metrics (genome
  fraction, NA50, misassemblies relative to the compare genome). When
  MUMmer's `dnadiff` is on PATH, `compare_dnadiff/out.report` adds a
  SNP/indel summary; both side reports are skipped silently when their
  binaries are missing. Works in two modes:
  *full pipeline* (`taco ... --compare X` runs steps 0–14 and the report
  lands in 14C), and *resume* (`taco ... -s 10,13,14 --compare X` after a
  prior full run still produces the report and a `compare` row in
  `final_result.csv`).
- **New `--final-fa <fasta>` flag.** Steps 13 and 14 now accept an
  externally produced final assembly via `--final-fa`. The flag is
  authoritative — it copies into `assemblies/final.merged.fasta` even
  when a stale copy is already present in the working directory. When
  `--final-fa` is omitted, TACO continues to fall back to
  `final_results/final.merged.fasta` (and `final_results/final_assembly.fasta`).
- **`--reference` is a full pipeline participant; `--compare` is fully
  passive.** ``EXCLUDED_FROM_REFINEMENT = {"compare"}`` controls every
  refinement-stage filter (auto-selector, fallback "first-available"
  picker, all-vs-all quickmerge candidate pool, telomere-pool input
  filter, polish/purge_dups loops). `--reference` therefore continues to
  contribute `*.telo.fasta` to the quickmerge candidate pool and is used
  as a chimera-detection alignment target — exactly as the original
  semantics imply. Only `--compare` is filtered out everywhere except QC
  and the new step 14C. (Use `--choose <assembler>` if you need to lock
  the backbone away from a high-scoring reference.)
- **Resume restoration covers steps 8–14.** Re-running individual steps
  after step-14 cleanup now restores normalized `assemblies/*.result.fasta`
  from `temp/assemblers/...`, the `assembly_info.csv` and merged metric
  CSVs from `final_results/`, telomere-pool FASTAs and provenance TSVs
  from `telomere_pool/`, and `final.merged.fasta` from
  `final_results/final.merged.fasta` (or from `--final-fa`).

### BUSCO: offline-first, configurable cache, and consistent across all steps

- **BUSCO trials are offline-first by default.** The previous default ran
  `--offline` only when the lineage was visible in a hard-coded location;
  otherwise it jumped straight to an online lookup that could timeout or
  attempt downloads on cluster nodes. Every BUSCO call (steps 8, 12, 13)
  now runs `--offline` first and only retries online when the offline
  attempt fails *and* the user has not passed `--busco-offline-only`.
- **New `--busco-download-path`.** Passed to BUSCO as `--download_path`
  so the cluster's central lineage cache is honored regardless of cwd.
  Falls back to the `BUSCO_DOWNLOAD_PATH` environment variable.
- **New `--busco-offline-only`.** Refuses any online lineage download.
  Replaces the older opt-in `STEP12_BUSCO_ALLOW_DOWNLOAD=1` flag (the env
  var still works, inverted: `STEP12_BUSCO_ALLOW_DOWNLOAD=0` disables the
  online fallback). Step 8's BUSCO loop now uses the same shared trial
  helper as steps 12/13, so all three sites benefit from the same logic.
- **`--taxon fungal` now defaults to `fungi_odb10`.** The previous default
  (`ascomycota_odb10`) is a sub-lineage; for a generic "fungal" choice
  TACO now picks the broad `fungi_odb10` lineage, matching the BUSCO
  recommendation and what the reference Venturiaceae paper used. Pass
  `--busco ascomycota_odb10` (or any other lineage) to override.

### Step 0 input QC: coverage estimator was always wrong

- **Bug fix: `step_00_input_qc` always reported the sample size as the
  total.** The estimator computed
  `bytes_per_read = fq_size / read_count; est_total_reads = fq_size / bytes_per_read`
  which is algebraically `read_count` — the per-read byte ratio used
  the *full* file size on both sides of the division. For a 12.57 GB ONT
  FASTQ the estimator reported ~100k reads / 459 Mb / 11.5× instead of
  the actual ~1.36M reads / 6.2 Gb / 155×. The estimator now (a) prefers
  `seqkit stats -T` for an exact count when `seqkit` is on PATH, and
  (b) falls back to a sampling estimator that tracks bytes consumed from
  the sample (or, for `.gz` inputs, the underlying compressed file's
  position). Verified on synthetic data: extrapolation is within ±0.5%
  on uncompressed and gzipped inputs.
- **Missing `gzip` import added** to `taco/steps.py` — the previous code
  referenced `gzip.open` without importing the module, which would have
  crashed on a `.fastq.gz` input.

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
- **purge_dups safety harmonized with genome-size budget.** The bp-drop
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
- **New** taxon-aware backbone scoring weights: fungi penalize BUSCO
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
completeness but did not penalize duplication.  An assembler with high BUSCO S
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
The README documents this dual-mode behavior.

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
