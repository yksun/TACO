"""Assembly purification: contaminant screening and chimera resolution.

This module holds the decision logic for TACO's purification step.  Every
function here is pure — it consumes per-contig statistics and parsed alignment
records and returns verdicts — so the science can be unit-tested without
minimap2, samtools, BUSCO, or a reference genome.  The orchestration that
produces the inputs lives in :mod:`taco.steps`.

Two problems are solved, and both were found on a real assembly of *Fusarium
tricinctum* MsR-QD66 (HiFi reads SRR33612568) checked against the published
Hi-C assembly GCA_050859235.1.

**Foreign sequence.**  purge_dups removes haplotypic duplication; it will not
remove another organism's genome.  A bacterial contig of 4.54 Mb at 19x against
the assembly's 313x modal depth, GC 66.6% against 47.7%, survived the whole
pipeline and inflated the published assembly length by 9% and its GC by 1.7
points.  Screening it out is easy; screening it out *without* also removing real
but atypical sequence is the hard part, and is what the guards below are for.

**Chimeric joins.**  When one assembler fuses two chromosomes through a
collapsed repeat, the result carries genuine telomere arrays at both outer ends
and looks unremarkable at the seam, so it passes a strict-T2T test and is
rewarded by the backbone score.  On the run above, Peregrine joined two
chromosomes through a collapsed ~7.9 kb tandem array at ~19x its assembled
depth; seven other assemblers and the Hi-C reference keep them apart, and *no
HiFi read spans the seam*.  Coverage across the seam does not dip — it spikes —
so a coverage-dip test sees nothing, and reads map inside the repeat rather than
clipping at its edge, so a clip-pileup test sees nothing either.  The decisive
evidence is the spanning-read count, corroborated by cross-assembler
disagreement.  This mirrors what Tigmint documents for repeat-mediated
mis-joins: coverage continuity is not evidence of a real join.

Design rules that follow from those two cases, and that the functions below
enforce:

1. **No single signal removes or breaks sequence.**  Established practice
   (BlobToolKit, NCBI FCS, Sanger ASCC, and the amosvalidate lineage of
   misassembly detectors) requires two independent evidence types before
   acting, because every individual axis is confounded.
2. **A robust z-score alone is not enough.**  On a genome whose real contigs
   cluster tightly, the robust sigma collapses and legitimate atypical
   chromosomes become outliers.  Every rule here requires a robust z-score
   *and* an absolute effect size.
3. **The baseline must survive its own contamination.**  A length-weighted
   median inverts host and contaminant once foreign sequence passes half the
   assembly, so the baseline is re-estimated on an iteratively trimmed core.
4. **Guards protect real-but-atypical sequence.**  Telomere-bearing contigs,
   high-depth organelles and collapsed repeats, and short contigs are never
   removed on composition alone.
"""

# ── robust statistics ────────────────────────────────────────────────────────

#: Scaling that makes the median absolute deviation a consistent estimator of
#: the standard deviation for normally distributed data.
MAD_TO_SIGMA = 1.4826

#: Robust z-score past which a value is considered an outlier.  This is the
#: Iglewicz-Hoaglin threshold, expressed in MAD-derived sigma units.
DEFAULT_Z_MIN = 3.5

#: Relative half-width of the depth cluster used to locate the core.  Real
#: chromosomes of one organism agree on depth to within a few percent; anything
#: needing more than a quarter of the modal value is a different population.
DEPTH_CLUSTER_TOL_REL = 0.25
#: Absolute half-width, in percentage points, of the GC cluster.
GC_CLUSTER_TOL_ABS = 3.0
#: Floors on the robust sigma, below which a z-score is not trustworthy.
#:
#: These are not cosmetic.  On the *F. tricinctum* assembly the real contigs
#: agree on depth to within about half a percent, so the unfloored MAD collapses
#: to ~1.5x on a 313x mode; a genuine chromosome 10% off the mode would then
#: score |z| = 20 and be condemned.  The GC MAD collapses to 0.66 points, which
#: is what makes the real accessory chromosome a 3.7-sigma outlier.  Flooring the
#: scale keeps the z-score meaningful on a tightly clustered assembly instead of
#: turning every small deviation into significance.
MIN_SIGMA_DEPTH_FRAC = 0.05   # as a fraction of the core depth
MIN_SIGMA_GC = 1.0            # percentage points


def median(values):
    """Plain median of *values*; 0.0 for an empty input."""
    s = sorted(float(v) for v in values)
    n = len(s)
    if n == 0:
        return 0.0
    if n % 2:
        return s[n // 2]
    return (s[n // 2 - 1] + s[n // 2]) / 2.0


def weighted_median(pairs):
    """Median of ``(value, weight)`` pairs, weighting by *weight*.

    Length-weighting stops a scatter of small aberrant contigs from moving the
    baseline, which is the common case.  It does *not* protect against a large
    contaminant, which is what :func:`robust_baseline` is for.
    """
    items = sorted((float(v), float(w)) for v, w in pairs
                   if w is not None and float(w) > 0)
    if not items:
        return 0.0
    half = sum(w for _, w in items) / 2.0
    run = 0.0
    for v, w in items:
        run += w
        if run >= half:
            return v
    return items[-1][0]


def mad_sigma(values, center):
    """MAD of *values* about *center*, scaled to a sigma estimate."""
    vals = [float(v) for v in values]
    if not vals:
        return 0.0
    return MAD_TO_SIGMA * median([abs(v - center) for v in vals])


def robust_z(value, center, sigma):
    """Signed deviation in robust sigma units; 0.0 when *sigma* is degenerate.

    Returning 0.0 rather than infinity for a zero sigma is deliberate: a
    degenerate spread means the data cannot support an outlier call, and the
    absolute-effect-size test is what should decide.
    """
    if sigma is None or sigma <= 0:
        return 0.0
    return (float(value) - float(center)) / float(sigma)


def robust_baseline(pairs, tol_abs=None, tol_rel=None):
    """Centre and robust sigma of the densest cluster, by weight.

    Args:
        pairs: ``(value, weight)`` per contig — weight is contig length.
        tol_abs: absolute half-width of a cluster, in the units of *value*.
        tol_rel: half-width as a fraction of the candidate centre.  Use this for
            depth, which is multiplicative and overdispersed; use *tol_abs* for
            composition.  When both are given the wider window wins.

    Returns:
        dict with ``center``, ``sigma``, ``n_used``, ``n_total`` and
        ``weight_frac`` — the share of total weight inside the winning cluster.

    A plain length-weighted median is the midpoint of the assembly by length, so
    once foreign sequence passes half the assembly the "baseline" becomes the
    contaminant and the host is measured against it.  Iterative trimming does
    not rescue this: the initial centre sits on the wrong cluster and the sigma
    is inflated by the very gap that should be resolved, so nothing is ever
    trimmed — textbook masking.

    Locating the heaviest cluster instead is both immune to that inversion and a
    better definition of "modal depth", which is what a coverage histogram peak
    means in purge_dups and in blobplot-style analyses.  ``weight_frac`` is
    reported so the caller can refuse to act when no cluster holds a clear
    majority, which is the signature of a genuine metagenome rather than a
    contaminated single-organism assembly.
    """
    items = [(float(v), float(w)) for v, w in pairs
             if v is not None and w is not None and float(w) > 0]
    if not items:
        return {"center": 0.0, "sigma": 0.0, "n_used": 0, "n_total": 0,
                "weight_frac": 0.0}

    total_w = sum(w for _, w in items)

    if tol_abs is None and tol_rel is None:
        # No cluster width requested: treat the whole distribution as the core,
        # which reduces this to a length-weighted median plus a robust sigma.
        members = items
    else:
        def window(candidate):
            widths = []
            if tol_abs is not None:
                widths.append(float(tol_abs))
            if tol_rel is not None:
                widths.append(abs(float(candidate)) * float(tol_rel))
            return max(widths)

        # Every observed value is a candidate centre; the winner is the one
        # whose window holds the most sequence.  Ties go to the heavier, then
        # higher, cluster so the result is deterministic.
        best = None
        for cand, _ in items:
            half = window(cand)
            members = [(v, w) for v, w in items if abs(v - cand) <= half]
            key = (sum(w for _, w in members), cand)
            if best is None or key > best[0]:
                best = (key, members)
        members = best[1] if best else items
    center = weighted_median(members)
    sigma = mad_sigma([v for v, _ in members], center)
    return {"center": center, "sigma": sigma, "n_used": len(members),
            "n_total": len(items),
            "weight_frac": (sum(w for _, w in members) / total_w) if total_w else 0.0}


def is_outlier(value, baseline, z_min, abs_min):
    """True when *value* is both a robust-z outlier and past an absolute floor.

    Returns ``(flag, z, deviation)``.  Requiring both conditions is what keeps a
    tightly clustered assembly from flagging its own atypical chromosomes: on
    the *F. tricinctum* run the real accessory chromosome scores z = -3.71 on GC,
    past the 3.5 threshold, and is retained only because its 2.44-point
    deviation is inside the 5-point floor.
    """
    if value is None:
        return False, 0.0, 0.0
    z = robust_z(value, baseline["center"], baseline["sigma"])
    dev = float(value) - float(baseline["center"])
    return (abs(z) >= z_min and abs(dev) >= abs_min), z, dev


# ── contaminant screening ────────────────────────────────────────────────────

#: A contig whose median depth falls below this fraction of the assembly's core
#: depth is anomalously shallow for a chromosome of the sequenced organism.
DEPTH_RATIO_MAX = 0.34
#: Depth above this multiple of core depth marks an organelle or a collapsed
#: repeat.  Both are real sequence, so this *protects* rather than accuses.
DEPTH_RATIO_ORGANELLE = 3.0
#: Absolute GC deviation, in percentage points, required alongside the z-score.
#: Intra-genomic GC within one eukaryote routinely spans tens of points, so a
#: small deviation is never sufficient however significant it looks.
GC_DEV_MIN = 5.0
#: Contigs shorter than this are too noisy in depth and composition to screen.
MIN_SCREEN_LEN = 10000
#: Composition alone may only condemn a contig at least this long.  Below it, a
#: GC outlier must be corroborated by depth.  An empirical null built by cutting
#: this assembly's real contigs into windows shows ~5% of legitimate 10 kb
#: windows deviate by more than 7.5 GC points.
MIN_COMPOSITION_LEN = 250000
#: Fraction of a contig at near-zero depth that counts as a corroborating
#: signal — foreign sequence recruits reads unevenly.
LOW_COV_FRAC_MIN = 0.10
#: Below this core weight fraction the assembly has no dominant organism and
#: automatic removal is refused.
MIN_CORE_WEIGHT_FRAC = 0.50
#: Automatic removal is refused outright if it would discard more than this
#: fraction of the assembly.  No screen should silently delete a quarter of a
#: genome; past this point the assembly needs a human, or ``--metagenome``.
MAX_REMOVAL_FRAC = 0.25


def screen_contaminants(contigs, depth_ratio_max=DEPTH_RATIO_MAX,
                        gc_dev_min=GC_DEV_MIN, min_len=MIN_SCREEN_LEN,
                        min_composition_len=MIN_COMPOSITION_LEN,
                        z_min=DEFAULT_Z_MIN, metagenome=False,
                        gc_null_p99=None):
    """Classify contigs as foreign, suspect, or clean.

    Args:
        contigs: iterable of dicts.  Required: ``contig``, ``length``.
            Optional: ``median_cov``, ``gc``, ``telomere_tier``, ``zero_bp``,
            ``low_bp``, ``aligned_to_compare`` (bool), ``busco_count`` (int).
            Extra keys are ignored.
        depth_ratio_max / gc_dev_min / min_len / min_composition_len / z_min:
            see the module constants.
        gc_null_p99: 99th percentile of GC deviation among windows of this
            assembly's own core sequence, when the caller has measured it.  It
            raises the absolute GC floor to what this genome legitimately does,
            which is the only defensible way to set that floor across genomes.
        metagenome: when True nothing is ever classified ``foreign``; every
            anomaly is reported as ``suspect``.  A metagenome has no single
            modal depth or composition, so the whole premise of the screen — one
            dominant organism — does not hold.

    Returns:
        ``(rows, baseline)``.  *rows* has one dict per contig with ``tier``
        (``foreign``, ``suspect``, ``organelle_candidate`` or ``clean``),
        the individual signals, the guards that fired, and a human-readable
        ``reason``.  *baseline* records what the comparison was made against,
        including ``core_weight_frac`` and ``removal_allowed``.

    Only ``foreign`` is ever removed by the caller, and only when at least two
    independent primary signals fire and no guard vetoes.  ``busco_count`` is
    recorded but deliberately never contributes to a removal decision: on the
    *F. tricinctum* run a real 3.08 Mb chromosome carries zero BUSCO genes while
    the bacterial contig carries three, so gene content does not separate them.
    """
    recs = [c for c in contigs if int(c.get("length") or 0) > 0]
    if not recs:
        return [], {"depth": None, "gc": None, "n_contigs": 0,
                    "core_weight_frac": 0.0, "removal_allowed": False,
                    "metagenome": bool(metagenome)}

    # A missing median_cov means "no coverage record", which is not 0x.  Such
    # contigs are excluded from the depth baseline and screened on composition
    # only, so absent data can never itself condemn a contig.
    depth_base = robust_baseline(
        [(c.get("median_cov"), c.get("length")) for c in recs
         if c.get("median_cov") is not None],
        tol_rel=DEPTH_CLUSTER_TOL_REL)
    gc_base = robust_baseline(
        [(c.get("gc"), c.get("length")) for c in recs
         if c.get("gc") is not None], tol_abs=GC_CLUSTER_TOL_ABS)

    # Floor the scales before any z-score is taken against them.
    depth_base["sigma"] = max(depth_base["sigma"],
                              MIN_SIGMA_DEPTH_FRAC * depth_base["center"])
    gc_base["sigma"] = max(gc_base["sigma"], MIN_SIGMA_GC)

    # A per-genome composition null, when the caller has measured one, raises the
    # absolute GC floor to whatever this assembly's own legitimate sequence does.
    # Cutting this genome's real contigs into 10 kb windows shows 1% of them
    # deviating by more than 24 GC points, so a fixed floor cannot be right for
    # every genome and every contig length.
    if gc_null_p99 is not None:
        gc_dev_min = max(gc_dev_min, float(gc_null_p99))

    core_frac = min(depth_base["weight_frac"] or 1.0,
                    gc_base["weight_frac"] or 1.0)
    removal_allowed = (not metagenome) and core_frac >= MIN_CORE_WEIGHT_FRAC

    rows = []
    for c in recs:
        rows.append(_classify_one(c, depth_base, gc_base, depth_ratio_max,
                                  gc_dev_min, min_len, min_composition_len,
                                  z_min, removal_allowed))

    # Final safety valve: never discard a large share of the assembly without a
    # human in the loop, however confident the per-contig evidence looks.
    total_bp = sum(r["length"] for r in rows)
    removed_bp = sum(r["length"] for r in rows if r["tier"] == "foreign")
    removal_frac = (removed_bp / float(total_bp)) if total_bp else 0.0
    if removal_allowed and removal_frac > MAX_REMOVAL_FRAC:
        removal_allowed = False
        for r in rows:
            if r["tier"] == "foreign":
                r["tier"] = "suspect"
                r["reason"] += (f"; removal refused: would discard "
                                f"{100 * removal_frac:.0f}% of the assembly, "
                                f"over the {100 * MAX_REMOVAL_FRAC:.0f}% cap — "
                                f"inspect manually or pass --metagenome")

    baseline = {"depth": depth_base, "gc": gc_base, "n_contigs": len(recs),
                "core_weight_frac": core_frac,
                "removal_allowed": removal_allowed,
                "removal_frac": round(removal_frac, 4),
                "metagenome": bool(metagenome)}

    rows.sort(key=lambda r: (r["tier"] != "foreign", -r["length"]))
    return rows, baseline


def _classify_one(c, depth_base, gc_base, depth_ratio_max, gc_dev_min,
                  min_len, min_composition_len, z_min, removal_allowed):
    name = c.get("contig")
    length = int(c["length"])
    depth = c.get("median_cov")
    gc = c.get("gc")
    tier_telo = (c.get("telomere_tier") or "").strip()

    depth = float(depth) if depth is not None else None
    gc = float(gc) if gc is not None else None

    depth_center = depth_base["center"]
    ratio = (depth / depth_center) if (depth is not None and depth_center > 0) else None

    depth_out, depth_z, _ = is_outlier(depth, depth_base, z_min, 0.0)
    gc_out, gc_z, gc_dev = is_outlier(gc, gc_base, z_min, gc_dev_min)

    # A depth outlier only accuses when it is *shallow*.  High depth means an
    # organelle or a collapsed repeat, both of which are real sequence.
    depth_low = bool(depth_out and ratio is not None and ratio < depth_ratio_max)
    depth_high = bool(ratio is not None and ratio >= DEPTH_RATIO_ORGANELLE)

    signals = []
    if depth_low:
        signals.append(f"depth {depth:.0f}x is {ratio:.2f} of core "
                       f"{depth_center:.0f}x (z={depth_z:.1f})")
    if gc_out:
        signals.append(f"GC {gc:.1f}% deviates {gc_dev:+.1f} points from core "
                       f"{gc_base['center']:.1f}% (z={gc_z:.1f})")

    # Corroborating signals: recorded, and able to substitute for a second
    # primary signal, but never sufficient alone.
    corroborating = []
    lowfrac = _low_cov_fraction(c, length)
    if lowfrac is not None and lowfrac >= LOW_COV_FRAC_MIN:
        corroborating.append(f"{100 * lowfrac:.0f}% of the contig at near-zero depth")
    if c.get("aligned_to_compare") is False:
        corroborating.append("no alignment to the compare genome")
    if c.get("busco_count") is not None:
        # Recorded for the report only.  See the docstring: gene content does
        # not separate foreign sequence from gene-poor real chromosomes.
        pass

    n_primary = int(depth_low) + int(gc_out)

    guards = []
    if length < min_len:
        guards.append(f"shorter than the {min_len:,} bp screening floor")
    if tier_telo == "strict_t2t":
        guards.append("telomere arrays at both ends")
    if depth_high:
        guards.append(f"depth {ratio:.1f}x core — organelle or collapsed repeat")
    if gc_out and not depth_low and length < min_composition_len:
        guards.append(f"composition-only call under {min_composition_len:,} bp")

    if depth_high and (gc_out or n_primary):
        tier = "organelle_candidate"
    elif guards:
        tier = "suspect" if (n_primary or corroborating) else "clean"
    elif n_primary >= 2 or (n_primary == 1 and corroborating):
        tier = "foreign" if removal_allowed else "suspect"
    elif n_primary == 1:
        tier = "suspect"
    else:
        tier = "clean"

    reason_bits = list(signals) + list(corroborating)
    if guards:
        reason_bits.append("retained: " + "; ".join(guards))
    if tier == "suspect" and not removal_allowed and n_primary >= 2:
        reason_bits.append("removal disabled (metagenome mode or no dominant "
                           "organism in the assembly)")

    return {
        "contig": name,
        "length": length,
        "median_cov": depth,
        "depth_ratio": round(ratio, 4) if ratio is not None else None,
        "depth_z": round(depth_z, 2),
        "gc": round(gc, 4) if gc is not None else None,
        "gc_deviation": round(gc_dev, 2) if gc is not None else None,
        "gc_z": round(gc_z, 2),
        "telomere_tier": tier_telo or "-",
        "busco_count": c.get("busco_count"),
        "n_primary_signals": n_primary,
        "n_corroborating": len(corroborating),
        "guards": "; ".join(guards) or "-",
        "tier": tier,
        "reason": "; ".join(reason_bits) or "no anomaly",
    }


def _low_cov_fraction(c, length):
    """Fraction of the contig recorded at zero or near-zero depth, or None."""
    zero = c.get("zero_bp")
    low = c.get("low_bp")
    if zero is None and low is None:
        return None
    try:
        worst = max(int(zero or 0), int(low or 0))
    except (TypeError, ValueError):
        return None
    if length <= 0:
        return None
    return min(1.0, worst / float(length))


def contaminant_removal_set(rows):
    """Names of contigs the screen classified ``foreign``."""
    return {r["contig"] for r in rows if r.get("tier") == "foreign"}


# ── chimera breakpoint estimation ────────────────────────────────────────────

#: An alignment must cover at least this many bases of the query to count as a
#: component of a candidate mis-join.
MIN_COMPONENT_BP = 100000
#: Two components must be this disjoint on the query — measured as overlap
#: relative to the shorter component — to imply a junction between them.
MAX_COMPONENT_OVERLAP_FRAC = 0.20


def consensus_breakpoint(voter_intervals, min_component_bp=MIN_COMPONENT_BP,
                         max_overlap_frac=MAX_COMPONENT_OVERLAP_FRAC):
    """Estimate a mis-join breakpoint from cross-assembler query intervals.

    Args:
        voter_intervals: ``{voter: {target_name: [(qstart, qend), ...]}}`` —
            where each target sequence of each voting assembly aligns on the
            suspect contig.  This is exactly what the concordance split vote
            already computes.
        min_component_bp: ignore targets covering less than this.
        max_overlap_frac: reject a voter whose two components overlap by more
            than this fraction of the shorter one, since overlapping components
            imply a repeat rather than a junction.

    Returns:
        dict with ``breakpoint`` (int or None), ``per_voter`` mapping voter to
        its own estimate, ``n_voters``, ``spread`` (max minus min estimate) and
        ``left_end`` / ``right_start`` — the tightest bracket any voter gives.

    Each voter that splits the contig into exactly two substantial, largely
    disjoint components contributes the midpoint of its interior gap; the
    consensus is the median of those.  On the *F. tricinctum* mis-join, four
    voters give 1,585,397 / 1,585,579 / 1,592,993 / 1,601,285 and the median
    lands 1,515 bp from the reference-derived seam midpoint, comfortably inside
    the 21 kb collapsed array where the exact base is arbitrary anyway.
    """
    per_voter = {}
    brackets = []
    for voter, targets in (voter_intervals or {}).items():
        spans = []
        for tname, ivs in (targets or {}).items():
            if not ivs:
                continue
            lo = min(s for s, _ in ivs)
            hi = max(e for _, e in ivs)
            covered = merged_span(ivs)
            if covered >= min_component_bp:
                spans.append((lo, hi, covered, tname))
        if len(spans) != 2:
            continue
        spans.sort()
        (lo1, hi1, cov1, _), (lo2, hi2, cov2, _) = spans
        overlap = max(0, hi1 - lo2)
        if overlap > max_overlap_frac * min(cov1, cov2):
            continue
        per_voter[voter] = int((hi1 + lo2) // 2)
        brackets.append((hi1, lo2))

    if not per_voter:
        return {"breakpoint": None, "per_voter": {}, "n_voters": 0,
                "spread": 0, "left_end": None, "right_start": None}

    ests = sorted(per_voter.values())
    return {
        "breakpoint": int(median(ests)),
        "per_voter": dict(per_voter),
        "n_voters": len(per_voter),
        "spread": ests[-1] - ests[0],
        # The tightest bracket is the most informative one: the largest
        # left-component end and the smallest right-component start.
        "left_end": max(b[0] for b in brackets),
        "right_start": min(b[1] for b in brackets),
    }


def merged_span(intervals):
    """Total length covered by *intervals* (list of half-open ``(start, end)``)."""
    if not intervals:
        return 0
    ordered = sorted(intervals)
    total = 0
    cur_s, cur_e = ordered[0]
    for s, e in ordered[1:]:
        if s <= cur_e:
            cur_e = max(cur_e, e)
        else:
            total += cur_e - cur_s
            cur_s, cur_e = s, e
    return total + (cur_e - cur_s)


# ── read-level confirmation of a candidate junction ──────────────────────────

#: Reads anchored on both flanks required to call a junction real.  DENTIST uses
#: three spanning reads with a unique anchor for the equivalent decision.
MIN_SPANNING_READS = 3
#: Anchor required on each side of the tested coordinate.
SPANNING_ANCHOR_BP = 1000
#: Spanning reads as a fraction of reads merely overlapping the position, below
#: which the junction is refuted even if the absolute count looks acceptable.
MIN_SPANNING_FRAC = 0.02


def spanning_read_verdict(n_spanning, n_overlapping, min_reads=MIN_SPANNING_READS,
                          min_frac=MIN_SPANNING_FRAC):
    """Decide whether reads support a candidate junction.

    Args:
        n_spanning: reads anchored on both sides of the breakpoint.
        n_overlapping: reads aligned across the position at all.  This is the
            local expectation the spanning count must be judged against.
        min_reads / min_frac: see the module constants.

    Returns:
        dict with ``verdict`` — ``"supported"``, ``"refuted"`` or
        ``"inconclusive"`` — plus ``spanning_frac`` and a ``reason``.

    The absolute count alone is not enough because depth varies by more than an
    order of magnitude across a collapsed repeat: at the *F. tricinctum* seam
    5,894 reads overlap the midpoint and *none* spans it, while a genuine locus
    on the same contig has 301 overlapping and 106 spanning.  Judging the
    spanning count as a fraction of the local pileup separates these cleanly.
    ``inconclusive`` is returned when there is too little data to say, and the
    caller must not break a contig on an inconclusive verdict.
    """
    n_spanning = int(n_spanning or 0)
    n_overlapping = int(n_overlapping or 0)

    if n_overlapping <= 0:
        return {"verdict": "inconclusive", "spanning_frac": None,
                "n_spanning": n_spanning, "n_overlapping": n_overlapping,
                "reason": "no reads aligned across the position"}

    frac = n_spanning / float(n_overlapping)

    if n_spanning >= min_reads and frac >= min_frac:
        return {"verdict": "supported", "spanning_frac": round(frac, 5),
                "n_spanning": n_spanning, "n_overlapping": n_overlapping,
                "reason": (f"{n_spanning} reads span the position "
                           f"({100 * frac:.1f}% of {n_overlapping} overlapping)")}

    # Refuting needs a pileup deep enough that "no spanning read" is meaningful.
    if n_overlapping >= max(min_reads * 10, 30):
        return {"verdict": "refuted", "spanning_frac": round(frac, 5),
                "n_spanning": n_spanning, "n_overlapping": n_overlapping,
                "reason": (f"only {n_spanning} of {n_overlapping} overlapping "
                           f"reads span the position ({100 * frac:.2f}%)")}

    return {"verdict": "inconclusive", "spanning_frac": round(frac, 5),
            "n_spanning": n_spanning, "n_overlapping": n_overlapping,
            "reason": (f"{n_overlapping} overlapping reads is too shallow to "
                       f"judge {n_spanning} spanning")}


# ── the chimera action decision ──────────────────────────────────────────────

def chimera_decision(contig, concordance, breakpoint_info, spanning,
                     action="split", protected_telomere=False):
    """Combine cross-assembler and read evidence into one action.

    Args:
        contig: contig name, for the reason string.
        concordance: verdict dict from :func:`taco.concordance.concordance_verdict`.
        breakpoint_info: result of :func:`consensus_breakpoint`.
        spanning: result of :func:`spanning_read_verdict`, or None when reads
            were unavailable.
        action: user-requested action — ``split``, ``replace``, ``report`` or
            ``off``.
        protected_telomere: whether this contig is in the protected T2T pool.

    Returns:
        dict with ``action_taken`` (``split``, ``replace``, ``report``, ``none``),
        ``breakpoint``, ``confidence`` and ``reason``.

    Sequence is only ever broken on *convergent* evidence: a cross-assembler
    majority that splits the contig, a breakpoint those assemblers agree on, and
    reads that fail to span it.  Any one of the three missing downgrades the
    outcome to a report, because breaking a real chromosome is worse than
    shipping a flagged one.  A contig carrying telomeres at both ends is not
    exempt — that is precisely the signature a fused pair of incomplete
    chromosome ends produces — but it does raise the bar to a read refutation
    rather than assembler disagreement alone.
    """
    verdict = (concordance or {}).get("verdict")
    n_split = int((concordance or {}).get("n_split") or 0)
    n_inf = int((concordance or {}).get("n_informative") or 0)
    bp = (breakpoint_info or {}).get("breakpoint")
    n_voters = int((breakpoint_info or {}).get("n_voters") or 0)
    span_verdict = (spanning or {}).get("verdict")

    def out(taken, confidence, reason):
        return {"contig": contig, "action_taken": taken, "breakpoint": bp,
                "confidence": confidence, "concordance_verdict": verdict,
                "n_split": n_split, "n_informative": n_inf,
                "breakpoint_voters": n_voters,
                "spanning_verdict": span_verdict or "not_tested",
                "reason": reason}

    if action == "off":
        return out("none", "not_assessed", "chimera handling disabled")

    if verdict != "mis_join_candidate":
        return out("none", "corroborated",
                   f"{n_split} of {n_inf} informative assemblies split this contig")

    if bp is None:
        return out("report", "low",
                   f"{n_split} of {n_inf} assemblies split this contig but they "
                   f"do not agree on a breakpoint")

    if action == "report":
        return out("report", "medium",
                   f"{n_split} of {n_inf} assemblies split this contig near "
                   f"{bp:,}; reporting only as requested")

    if span_verdict == "supported":
        return out("none", "high",
                   f"{n_split} of {n_inf} assemblies split this contig but reads "
                   f"span {bp:,} — the join is read-supported and is kept "
                   f"({(spanning or {}).get('reason', '')})")

    if span_verdict == "refuted":
        return out(action, "high",
                   f"{n_split} of {n_inf} assemblies split this contig near "
                   f"{bp:,} and reads do not span it "
                   f"({(spanning or {}).get('reason', '')})")

    # No read evidence either way.  Assembler disagreement alone is enough to
    # act only when it is unanimous and the breakpoint is well determined.
    unanimous = n_inf >= 3 and n_split == n_inf and n_voters >= 2
    if unanimous and not protected_telomere:
        return out(action, "medium",
                   f"all {n_inf} informative assemblies split this contig near "
                   f"{bp:,}; no read evidence available")

    return out("report", "low",
               f"{n_split} of {n_inf} assemblies split this contig near {bp:,}; "
               f"read evidence unavailable"
               + (" and the contig is telomere-protected" if protected_telomere else "")
               + " — reporting without breaking it")


def apply_split(name, seq, breakpoint, suffix_a="a", suffix_b="b"):
    """Split one contig into two at *breakpoint*.

    Returns a list of ``(name, sequence)``, or the original single record when
    the breakpoint is not strictly interior.  Naming keeps the parent visible so
    provenance stays traceable in the FASTA itself.
    """
    if breakpoint is None or breakpoint <= 0 or breakpoint >= len(seq):
        return [(name, seq)]
    return [(f"{name}_{suffix_a}", seq[:breakpoint]),
            (f"{name}_{suffix_b}", seq[breakpoint:])]
