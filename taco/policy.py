"""Assembly representation policy: how an assembly is scored and what may be removed.

TACO can pursue two different objectives, and they are not the same objective
with different thresholds — they disagree about what a good assembly *is*.

``primary``
    One nonredundant representative sequence per chromosome.  A duplicated
    single-copy ortholog is evidence of retained redundancy, so BUSCO
    duplication is a real negative signal and collapsing steps are wanted.

``full``
    Full sequence representation: divergent alternate/haplotype sequence is
    kept rather than collapsed.  In a heterozygous, diploid, or dikaryotic
    sample a high BUSCO duplicated percentage is the expected consequence of
    representing both copies, not an error, so it must not be penalised as one.

``full`` is **not** a phased assembly.  TACO does not order or orient contigs
into haplotypes and does not validate haplotype consistency, so nothing here may
be described as phased, haplotype-resolved, or hap1/hap2.  The mode retains
alternate sequence; it does not assign it.

Everything in this module is a pure function of metrics and configuration, so the
policy can be unit-tested without an assembly, a reference, or any external tool.
Decisions that consume it live in :mod:`taco.steps`.
"""

# ── Scoring profiles ─────────────────────────────────────────────────────────
#
# Weights are deliberately kept on the same numeric scale as the pre-v1.5.0
# formula so that scores remain comparable across releases and so the primary
# profile reproduces historical selections.  The taxon-specific overrides that
# existed before are applied on top of the profile, not replaced by it.

#: Which BUSCO quantity the profile rewards.  ``single`` is the nonredundant
#: objective; ``complete`` (S + D) is the full-representation objective, under
#: which a locus present twice is still a locus present.
BUSCO_TERM_SINGLE = "single"
BUSCO_TERM_COMPLETE = "complete"

#: Weights of the pre-v1.5.0 formula, kept as a documented record and as the
#: oracle for the backward-compatibility test.  Not selectable: v1.5.0
#: deliberately changed how a backbone is chosen (see REFINEMENT-AWARE below).
LEGACY_PRIMARY_WEIGHTS = {
    "busco_term": BUSCO_TERM_SINGLE,
    "w_busco": 1000, "w_busco_d": 500, "w_t2t": 300, "w_single_tel": 150,
    "w_merqury_comp": 200, "w_merqury_qv": 20, "w_contigs": 30,
    "w_contig_density": 0, "w_n50": 150, "w_size_penalty": 500,
}

# ── REFINEMENT-AWARE SELECTION ───────────────────────────────────────────────
#
# The backbone is chosen BEFORE refinement runs, so the score should weigh what
# refinement cannot change more heavily than what it demonstrably fixes.
# Measured on the real Puccinia triticina run, backbone (ipa) -> delivered:
#
#     BUSCO C (%)            96.3  -> 96.2     UNFIXABLE, only degraded
#     Merqury completeness   91.74 -> 91.51    UNFIXABLE, only degraded
#     N50 (Mb)                2.28 ->  2.48    barely movable
#     Merqury QV             40.50 -> 47.85    fixable: polishing
#     BUSCO D (%)              6.9 ->   5.3    fixable: purge_dups
#     strict T2T                 1 ->     4    fixable: telomere-pool rescue
#     contigs                  329 ->   117    fixable downward
#
# The pre-v1.5.0 weights had this backwards: BUSCO_D was charged 600/point --
# the heaviest term after BUSCO itself -- for duplication that step 12H removes
# anyway, while contiguity got 30/contig and consensus accuracy 20/QV-point.
# It therefore rejected a 99-contig / N50 4.08 Mb / QV 71.7 candidate in favour
# of a 329-contig / N50 2.28 Mb / QV 40.5 one, and delivered an assembly less
# contiguous and 24 QV points worse than the candidate it passed over.
#
# MERQURY COMPLETENESS MEANS OPPOSITE THINGS IN THE TWO MODES, which is why it
# is not simply "a quality metric".  It is the fraction of read k-mers present
# in the assembly.  A collapsed haploid representation of a heterozygous sample
# CANNOT contain both haplotypes' k-mers, so:
#   * in full mode it measures the objective directly and dominates the score;
#   * in primary mode a value above the haploid ceiling is evidence of retained
#     redundancy rather than of quality, so it is weighted lightly and BUSCO C
#     -- which a haploid representation can legitimately reach 100% of -- is the
#     content signal instead.
# On the real data the collapsed candidates sit at 91.7-93.8% completeness and
# the haplotype-retaining ones at 98.7-99.2%, so this distinction decides which
# assembly each mode selects.

SCORING_PROFILES = {
    "primary": {
        # BUSCO C, not S: S is partly a function of duplication, which
        # purge_dups removes, so rewarding S pre-purge rewards work already
        # done and charging D pre-purge charges it a second time.
        "busco_term": BUSCO_TERM_COMPLETE,
        "w_busco": 1000,
        # Duplication is the primary objective's defect but it is FIXABLE, so it
        # informs the choice without dominating it.  Taxon overrides scale it.
        "w_busco_d": 150,
        # Fixable by telomere-pool rescue (1 -> 4 T2T observed), so modest.
        "w_t2t": 120,
        "w_single_tel": 50,
        # Light: see the note above -- excess k-mer completeness is redundancy
        # under this objective, not quality.
        "w_merqury_comp": 150,
        # Partly fixable by polishing (+7.35 observed) but a 31-point QV spread
        # is a ~1300x difference in consensus error, which 20/point ignored.
        "w_merqury_qv": 120,
        "w_contigs": 60,          # unfixable upward; taxon overrides scale it
        "w_contig_density": 0,
        "w_n50": 600,             # unfixable; applied to log10(N50)
        "w_size_penalty": 500,
        "ploidy_expectation": 1.0,
        "description": "nonredundant primary representation",
    },
    "full": {
        "busco_term": BUSCO_TERM_COMPLETE,
        "w_busco": 1000,
        # Duplication is approximately neutral: not rewarded, not punished.
        # In a full representation a duplicated ortholog is the expected
        # consequence of the objective.
        "w_busco_d": 0,
        "w_t2t": 120,
        "w_single_tel": 50,
        # DOMINANT: k-mer completeness is the direct measurement of "did the
        # assembly keep all the read-supported sequence", which IS this mode's
        # objective.  Without it the mode has no positive signal for haplotype
        # retention at all -- BUSCO C saturates near 96% for every serious
        # candidate, and duplication and size are deliberately not scored -- and
        # it selects a collapsed 147 Mb assembly over a genuine 254 Mb one.
        # The value is not knife-edge: on the real data the same assembly wins
        # anywhere from 1500 to 4000, so this only widens the margin over the
        # best collapsed candidate (0.5% at 1500, 2.1% here) rather than
        # deciding the answer.
        "w_merqury_comp": 2500,
        "w_merqury_qv": 120,
        "w_contigs": 0,               # absolute count not charged here
        # Charged at the same per-contig rate as primary, but against a
        # haploid-sized assembly, so an assembly that legitimately carries 2n
        # contigs is not charged twice for doing so.
        "w_contig_density": 30,
        "w_n50": 600,
        # Size does not RANK in full mode.  No assembler declares an expected
        # total size for a heterozygous sample: hifiasm and Verkko take no size
        # at all, Flye's is optional, and Canu's genomeSize is a coverage knob
        # (corOutCoverage, mhapSensitivity, NG50 logging).  Where a declared
        # size IS compared against a finished assembly, NCBI uses a wide
        # taxon-derived band anchored on the haploid value to FLAG for human
        # review, never to rank.  The band is reporting and advisory only.
        "w_size_penalty": 0,
        "ploidy_expectation": None,
        "description": "full sequence representation, alternate sequence retained",
    },
}

#: Taxon overrides, applied to whichever profile is active.  These predate the
#: profiles and are preserved so primary-mode selections do not move.
TAXON_OVERRIDES = {
    # w_busco_d keeps the RELATIVE ordering the taxon presets always had
    # (fungal strictest, plant most permissive because polyploidy inflates D),
    # rescaled to the refinement-aware level: duplication is fixable, so it no
    # longer outweighs contiguity.  w_contigs keeps its preset ordering too.
    "fungal":     {"w_busco_d": 180, "w_t2t": 140, "w_n50": 600, "w_contigs": 60},
    "plant":      {"w_busco_d": 90,  "w_t2t": 80,  "w_n50": 600, "w_contigs": 100},
    "vertebrate": {"w_busco_d": 150, "w_t2t": 80,  "w_n50": 800, "w_contigs": 80},
    "animal":     {"w_busco_d": 150, "w_t2t": 80,  "w_n50": 800, "w_contigs": 80},
    "insect":     {"w_busco_d": 150, "w_t2t": 120, "w_n50": 600, "w_contigs": 60},
}

#: Modes that name a scoring profile and therefore produce one deliverable.
DELIVERABLE_MODES = ("primary", "full")
#: ``both`` is a DRIVER mode, not a scoring profile: it runs the deliverable
#: steps once per profile and emits both genomes from one set of assemblies.
MODE_BOTH = "both"

#: ``both`` is the default: the two objectives disagree about what a good
#: assembly is, and which one a given project wants is not knowable from the
#: reads, so TACO emits both rather than silently choosing for the user.
DEFAULT_MODE = MODE_BOTH
VALID_MODES = DELIVERABLE_MODES + (MODE_BOTH,)


def resolve_weights(mode=DEFAULT_MODE, taxon="other"):
    """Profile weights for *mode*, with *taxon* overrides applied.

    A duplication penalty is only meaningful where duplication is unwanted, so a
    taxon override for ``w_busco_d`` is ignored in a profile that has set it to
    zero on purpose.  Without that guard, ``--taxon fungal`` would silently
    reimpose a 600-point duplication penalty on full mode and defeat it.
    """
    if mode == MODE_BOTH:
        raise ValueError(
            "'both' names no single weight set; it produces one deliverable per "
            "profile. Resolve weights per deliverable mode instead "
            f"({', '.join(DELIVERABLE_MODES)}).")
    if mode not in SCORING_PROFILES:
        raise ValueError(f"unknown assembly mode {mode!r}; expected one of "
                         f"{', '.join(VALID_MODES)}")
    w = dict(SCORING_PROFILES[mode])
    for key, val in (TAXON_OVERRIDES.get(taxon) or {}).items():
        if key == "w_busco_d" and w.get("w_busco_d", 0) == 0:
            continue
        if key == "w_contigs" and w.get("w_contig_density", 0) > 0:
            # This profile charges contig density instead; a taxon override for
            # the absolute count would silently restore the bias against a
            # haplotype-retaining assembly.
            continue
        w[key] = val
    return w


# ── Size expectation ─────────────────────────────────────────────────────────
#
# The size term is the part of the score that primary and full modes cannot
# share, and getting it wrong is the main way full mode fails.
#
# Assemblers take -g as the HAPLOID size (Canu genomeSize, Flye --genome-size),
# so in primary mode the delivered assembly should approximate it and a
# point-deviation penalty is correct.  In full mode a heterozygous or dikaryotic
# assembly legitimately approaches a multiple of the haploid size, and the same
# penalty then punishes the objective: at -g 127m, a 254 Mb dikaryon takes a
# ~500-point hit purely for being what the user asked for.
#
# We therefore use a TOLERANCE BAND rather than a point in full mode: anything
# between the haploid size and the ploidy expectation is unpenalised, and only
# distance OUTSIDE the band is charged.  The band is deliberately not a bare
# 2x multiplication — the multiple is a parameter, because real dikaryons are
# not exactly 2x and the user may know their ploidy.
#
# In full mode the band does NOT contribute to the score (w_size_penalty = 0);
# it exists to populate the report and the advisory.  The guard against an
# inflated assembly is structural — contig density and N50 — because a 287 Mb /
# 46,812-contig artifact sits INSIDE any band wide enough to admit a real
# 254 Mb dikaryon, so size cannot be the discriminator.

#: Upper multiple of the haploid size treated as legitimate in full mode when
#: the user has not said otherwise.  2.0 covers a diploid or dikaryon; the
#: slack above it absorbs the fact that real dikaryons are not exactly 2x.
DEFAULT_FULL_PLOIDY = 2.0
#: Fractional slack added above the ploidy multiple before any penalty applies.
PLOIDY_BAND_SLACK = 0.15


def size_expectation(total_len, expected_haploid, mode=DEFAULT_MODE,
                     ploidy=None, band_slack=PLOIDY_BAND_SLACK):
    """Relative size deviation to be charged, and the band it was judged against.

    Args:
        total_len: assembly length in bp.
        expected_haploid: the user's ``-g`` value in bp (haploid).
        mode: ``primary`` or ``full``.
        ploidy: upper multiple of the haploid size regarded as legitimate in
            full mode.  ``None`` uses :data:`DEFAULT_FULL_PLOIDY`.
        band_slack: fractional slack above the upper bound.

    Returns:
        ``(deviation, low, high)`` — *deviation* is the fraction to charge
        (0.0 when the assembly is inside the band), and *low*/*high* describe
        the band for reporting.  Primary mode returns a degenerate band where
        low == high == expected_haploid, reproducing the historical behaviour
        exactly.
    """
    if not expected_haploid or expected_haploid <= 0 or not total_len or total_len <= 0:
        return 0.0, None, None

    if mode != "full":
        return abs(total_len - expected_haploid) / expected_haploid, \
            expected_haploid, expected_haploid

    p = DEFAULT_FULL_PLOIDY if ploidy is None else float(ploidy)
    low = expected_haploid
    high = expected_haploid * p * (1.0 + band_slack)
    if low <= total_len <= high:
        return 0.0, low, high
    if total_len < low:
        return (low - total_len) / low, low, high
    return (total_len - high) / high, low, high


def _contig_density(contigs, total_len):
    """Contigs per Mb of assembled sequence, or None if the length is unknown."""
    if not total_len or total_len <= 0:
        return None
    return contigs / (total_len / 1e6)


def _contig_density_penalty(contigs, total_len, expected_haploid, w):
    """Fragmentation charge, normalised to haploid size.

    Primary mode charges the absolute contig count.  That is correct when the
    target is one collapsed copy of the genome, but it structurally penalises a
    haplotype-retaining assembly, which carries roughly ploidy-times as many
    contigs for exactly the reason full mode exists.  So full mode charges
    contig DENSITY at the same per-contig rate, evaluated against a
    haploid-sized assembly:

        density (contigs/Mb) x rate (per contig) x haploid size (Mb)

    which is the cost the absolute penalty would have levied on an assembly of
    haploid length at the same density.  Fragmentation is therefore priced
    identically in both modes; only the double-charge for ploidy is removed.

    With no declared genome size the normaliser falls back to the assembly's own
    length, which reduces the expression exactly to the absolute count penalty.
    """
    rate = w.get("w_contig_density", 0)
    if not rate:
        return 0.0
    density = _contig_density(contigs, total_len)
    if density is None:
        return 0.0
    norm_mb = (expected_haploid / 1e6) if expected_haploid else (total_len / 1e6)
    return density * rate * norm_mb


#: Smallest span treated as a plausible chromosome when bounding the telomere
#: reward.  Deliberately far below any real chromosome in the genomes TACO
#: targets (rust chromosomes average ~7 Mb), so the bound below only ever fires
#: on counts that are impossible by a wide margin rather than merely surprising.
MIN_PLAUSIBLE_CHROMOSOME_BP = 1_000_000


def telomere_reward_cap(expected_haploid, mode=DEFAULT_MODE, ploidy=None):
    """Largest number of telomere-bearing contigs worth rewarding, or None.

    A contig has two ends, so it can carry at most two telomeres, and a genome
    holds 2 x n_chromosomes x ploidy telomeres in total.  The count of contigs
    that can legitimately carry one is therefore bounded by the chromosome
    count -- which TACO does not know, so it is bounded from above here by
    assuming chromosomes no smaller than :data:`MIN_PLAUSIBLE_CHROMOSOME_BP`.

    Without this bound the reward is unbounded in a quantity that biology
    bounds, and a shattered assembly accumulates repeat-derived telomere-like
    ends without limit.  On real Puccinia triticina data one candidate reported
    1,668 single-end-strong telomere contigs in an 18-chromosome genome (at most
    36 are possible, 72 as a dikaryon), earning three times its entire BUSCO
    term; only its contig count kept it from winning.

    Returns None when no genome size was declared, in which case no bound can be
    derived and the historical unbounded behaviour is kept.
    """
    if not expected_haploid or expected_haploid <= 0:
        return None
    p = 1.0 if mode != "full" else (DEFAULT_FULL_PLOIDY if ploidy is None
                                    else float(ploidy))
    n_chromosomes = (expected_haploid * p) / MIN_PLAUSIBLE_CHROMOSOME_BP
    # x2: a chromosome broken into two pieces leaves each piece one telomere.
    return max(1.0, 2.0 * n_chromosomes)


def partially_available_metrics(all_metrics, optional=("merqury_qv", "merqury_comp")):
    """Optional metrics present for SOME candidates but not all.

    A metric that some candidates lack is not comparable across the set.
    Scoring the absentees as zero does not make the term neutral -- zero is the
    minimum, so it is the maximum penalty, and raising a weight makes the
    distortion proportionally worse.  On a real run two of nine candidates had
    no Merqury QV, which cost them 6,600 points against a peer at QV 55.

    Returns the names of such metrics so the caller can drop those terms for
    every candidate rather than penalising the ones with missing data.
    """
    def _numeric(v):
        # A non-finite or non-numeric value (Merqury reports "+inf" when it
        # finds no unsupported k-mers) is not scoreable, so it counts as absent.
        try:
            f = float(v)
        except (TypeError, ValueError):
            return False
        return f == f and abs(f) != float("inf") and f != 0

    partial = set()
    for key in optional:
        seen = [m.get(key) for m in all_metrics]
        have = [v for v in seen if _numeric(v)]
        if have and len(have) < len(seen):
            partial.add(key)
    return partial


def score_assembly(metrics, mode=DEFAULT_MODE, taxon="other",
                   expected_haploid=0, ploidy=None, unavailable=()):
    """Composite selection score for one candidate assembly.

    Args:
        metrics: mapping with any of ``busco_s``, ``busco_d``, ``busco_c``,
            ``t2t``, ``single_tel``, ``merqury_comp``, ``merqury_qv``,
            ``contigs``, ``n50``, ``total_len``.  Missing values count as zero,
            which is how an unavailable metric (e.g. Merqury not installed)
            contributes nothing rather than penalising the candidate.
        mode / taxon: select the weight set.
        expected_haploid / ploidy: passed to :func:`size_expectation`.

    Returns:
        ``(score, contributions)``.  *contributions* maps each named component
        to its signed point value, so a report can state why an assembly won
        instead of only that it did.
    """
    import math

    w = dict(resolve_weights(mode, taxon))
    # A metric that is not comparable across the candidate set is dropped for
    # everyone; see partially_available_metrics.
    for key in unavailable:
        if key == "merqury_qv":
            w["w_merqury_qv"] = 0
        elif key == "merqury_comp":
            w["w_merqury_comp"] = 0
    def g(k):
        try:
            v = float(metrics.get(k) or 0)
        except (TypeError, ValueError):
            return 0.0                 # e.g. Merqury's "+inf (no unsupported k-mers)"
        return 0.0 if (v != v or abs(v) == float("inf")) else v

    busco_s, busco_d, contigs = g("busco_s"), g("busco_d"), g("contigs")
    n50, total_len = g("n50"), g("total_len")
    # BUSCO C is preferred when reported; S + D reconstructs it otherwise, which
    # is what makes full mode work on tables that only carry S and D.
    busco_c = g("busco_c") or (busco_s + busco_d)
    busco_used = busco_c if w["busco_term"] == BUSCO_TERM_COMPLETE else busco_s

    deviation, low, high = size_expectation(
        total_len, expected_haploid, mode=mode, ploidy=ploidy)

    # Telomere counts are rewarded, and biology bounds them; see
    # telomere_reward_cap.  Counts past the bound are false positives, so they
    # are not rewarded further rather than being penalised.
    telo_cap = telomere_reward_cap(expected_haploid, mode=mode, ploidy=ploidy)
    t2t_raw, single_raw = g("t2t"), g("single_tel")
    if telo_cap is None:
        t2t_used, single_used = t2t_raw, single_raw
    else:
        t2t_used = min(t2t_raw, telo_cap)
        single_used = min(single_raw, telo_cap)

    contributions = {
        "busco": busco_used * w["w_busco"],
        "duplication_penalty": -busco_d * w["w_busco_d"],
        "t2t": t2t_used * w["w_t2t"],
        "single_telomere": single_used * w["w_single_tel"],
        "merqury_completeness": g("merqury_comp") * w["w_merqury_comp"],
        "merqury_qv": g("merqury_qv") * w["w_merqury_qv"],
        "contig_count_penalty": -contigs * w.get("w_contigs", 0),
        "contig_density_penalty": -_contig_density_penalty(
            contigs, total_len, expected_haploid, w),
        "contiguity": (math.log10(n50) * w["w_n50"]) if n50 > 0 else 0.0,
        "size_penalty": -deviation * w["w_size_penalty"],
    }
    meta = {"busco_term": w["busco_term"], "size_band_low": low,
            "size_band_high": high, "size_deviation": deviation,
            "size_ranks": bool(w.get("w_size_penalty", 0)),
            "contigs_per_mb": _contig_density(contigs, total_len),
            "telomere_cap": telo_cap,
            "telomere_capped": bool(telo_cap is not None
                                    and (t2t_raw > telo_cap or single_raw > telo_cap))}
    return sum(contributions.values()), {**contributions, "_meta": meta}


# ── What may be removed ──────────────────────────────────────────────────────

class AssemblyPolicy:
    """Which sequence-removing operations are permitted under a mode.

    Only steps that remove sequence *because it looks redundant* are gated.
    Contamination screening and chimera resolution are never gated: foreign
    sequence and unsupported joins are errors under either objective, and they
    are identified from independent evidence (composition, read depth, read
    spanning) rather than from resembling a second haplotype.

    ``both`` is accepted here so a run can carry one policy object, but it
    answers no removal question itself: the two objectives disagree about what
    may be removed, so the caller must ask a sub-policy.  Use
    :meth:`sub_policies`.
    """

    def __init__(self, mode=DEFAULT_MODE, taxon="other", user_no_purge_dups=False,
                 ploidy=None):
        if mode not in VALID_MODES:
            raise ValueError(f"unknown assembly mode {mode!r}; expected one of "
                             f"{', '.join(VALID_MODES)}")
        self.mode = mode
        self.taxon = taxon
        self.ploidy = ploidy
        self.user_no_purge_dups = bool(user_no_purge_dups)

    @property
    def is_both(self):
        return self.mode == MODE_BOTH

    @property
    def is_full(self):
        return self.mode == "full"

    def sub_policies(self):
        """One policy per deliverable this run will produce.

        For a single-mode run that is just this policy; for ``both`` it is one
        per profile, in the order the deliverables are produced.
        """
        if not self.is_both:
            return [self]
        return [AssemblyPolicy(mode=m, taxon=self.taxon,
                              user_no_purge_dups=self.user_no_purge_dups,
                              ploidy=self.ploidy)
                for m in DELIVERABLE_MODES]

    def _reject_both(self, what):
        if self.is_both:
            raise ValueError(
                f"{what} is not defined for mode 'both': the primary and full "
                f"objectives disagree about it. Ask a sub_policies() member.")

    @property
    def purge_dups_enabled(self):
        """purge_dups removes haplotigs by design, so full mode declines it.

        An explicit ``--no-purge-dups`` always wins; it can disable the step but
        never re-enable it against the mode.
        """
        self._reject_both("purge_dups")
        if self.user_no_purge_dups:
            return False
        return not self.is_full

    @property
    def collapse_redundancy_enabled(self):
        """Self-dedup and containment filtering also delete alternate copies.

        Two near-identical contigs are exactly what a retained haplotype pair
        looks like, so these passes are skipped in full mode.
        """
        self._reject_both("redundancy collapsing")
        return not self.is_full

    @property
    def contamination_screen_enabled(self):
        return True          # never gated by representation mode

    @property
    def chimera_resolution_enabled(self):
        return True          # never gated by representation mode

    def describe(self):
        if self.is_both:
            parts = "; ".join(
                f"{m} ({SCORING_PROFILES[m]['description']})"
                for m in DELIVERABLE_MODES)
            return f"Assembly representation mode: {MODE_BOTH} — {parts}"
        prof = SCORING_PROFILES[self.mode]
        return (f"Assembly representation mode: {self.mode} "
                f"({prof['description']})")


# ── Advisory ─────────────────────────────────────────────────────────────────

#: Thresholds for the advisory below.  These describe a pattern worth telling
#: the user about; they are not a ploidy call and must never be reported as one.
ADVISORY_MIN_BUSCO_D = 20.0
ADVISORY_MIN_MERQURY_COMP = 95.0
ADVISORY_MIN_SIZE_RATIO = 1.4


def haplotype_retention_advisory(metrics, expected_haploid, mode=DEFAULT_MODE):
    """Advise when metrics look like retained alternate sequence, without concluding it.

    High k-mer completeness together with strong BUSCO duplication and a
    substantially larger assembly is the signature of a representation that has
    kept both copies.  It is also, less often, the signature of artificial
    duplication — which is why this returns advice and not a diagnosis.  Nothing
    here infers ploidy, and the text must not claim the sample is diploid or
    dikaryotic.
    """
    if mode in ("full", MODE_BOTH):
        # Either the user asked for a full representation, or 'both' is already
        # delivering one alongside the primary; advising it would be noise.
        return None
    busco_d = float(metrics.get("busco_d") or 0)
    comp = float(metrics.get("merqury_comp") or 0)
    total_len = float(metrics.get("total_len") or 0)
    if not (expected_haploid and total_len):
        return None
    ratio = total_len / float(expected_haploid)
    if not (busco_d >= ADVISORY_MIN_BUSCO_D
            and comp >= ADVISORY_MIN_MERQURY_COMP
            and ratio >= ADVISORY_MIN_SIZE_RATIO):
        return None
    return (
        "High k-mer completeness coincides with strong BUSCO duplication and "
        f"increased assembly size ({comp:.1f}% complete, {busco_d:.1f}% "
        f"duplicated, {ratio:.2f}x the expected haploid size).\n"
        "This may indicate retention of alternate haplotypes or duplicated "
        "sequence rather than a superior primary assembly. These metrics alone "
        "do not establish the ploidy of the sample.\n"
        "Consider --assembly-mode full if a full sequence representation is "
        "desired.")
