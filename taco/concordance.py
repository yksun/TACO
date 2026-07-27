"""Cross-assembler concordance, mis-join detection, and contaminant screening.

A multi-assembler pipeline has evidence no single-assembly tool has: when one
assembler joins two sequences that every other assembler keeps apart, that join
is a chimera candidate.  This matters specifically for telomere-based scoring,
because fusing two *incomplete* chromosome ends produces a contig with genuine
telomeres at both outer ends and nothing anomalous at the seam, so it passes a
strict-T2T test and is rewarded by the backbone score.  Such a join cannot be
detected by inspecting the contig alone; it is only visible as disagreement
between assemblers.

The functions here are deliberately pure: they consume parsed alignment records
and per-contig statistics and return verdicts, so the decision logic can be
unit-tested without minimap2, BUSCO, or a reference genome.  The orchestration
that produces their inputs lives in ``taco.steps``.

Three independent checks are provided:

``concordance_verdict``
    Does another assembly split this contig?  Majority vote across assemblers.
``detect_fusions``
    In ``--compare`` mode, do two reference sequences land on one contig?
``screen_contaminants``
    Is a contig anomalous in read depth and/or GC relative to the assembly?
"""

# ── Cross-assembler split voting ─────────────────────────────────────────────

#: A single target must cover at least this fraction of the query contig for
#: the query to be considered reproduced intact by that assembly.
DOMINANT_FRAC = 0.90
#: Each of two or more targets must cover at least this fraction of the query
#: for the query to be considered split by that assembly.
SEGMENT_FRAC = 0.20
#: Total aligned coverage below this fraction leaves the assembly with no
#: opinion (too diverged, too fragmented, or genuinely missing the region).
INFORMATIVE_FRAC = 0.60
#: Minimum number of assemblies that must vote "split" before a contig is
#: called a mis-join candidate, regardless of the majority.
MIN_SPLIT_VOTES = 2


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


def split_vote(query_len, target_intervals,
               dominant_frac=DOMINANT_FRAC, segment_frac=SEGMENT_FRAC,
               informative_frac=INFORMATIVE_FRAC):
    """Decide whether one assembly reproduces a contig intact or splits it.

    Args:
        query_len: length of the contig being tested.
        target_intervals: ``{target_name: [(qstart, qend), ...]}`` — where each
            target sequence of the other assembly aligns *on the query*.
        dominant_frac / segment_frac / informative_frac: see module constants.

    Returns:
        One of ``"intact"``, ``"split"``, or ``"uninformative"``.

    A contig is ``split`` when two or more targets each cover a substantial and
    largely distinct portion of it while no single target covers most of it.
    Requiring both conditions avoids calling a split on the basis of a repeat
    that happens to align in several places.
    """
    if query_len <= 0 or not target_intervals:
        return "uninformative"

    covers = {t: merged_span(iv) for t, iv in target_intervals.items()}
    total = merged_span([iv for ivs in target_intervals.values() for iv in ivs])
    if total < informative_frac * query_len:
        return "uninformative"

    best = max(covers.values()) if covers else 0
    if best >= dominant_frac * query_len:
        return "intact"

    substantial = [t for t, c in covers.items() if c >= segment_frac * query_len]
    if len(substantial) >= 2:
        return "split"
    return "uninformative"


def concordance_verdict(votes, min_split_votes=MIN_SPLIT_VOTES):
    """Combine per-assembly votes for one contig into a verdict.

    Args:
        votes: ``{assembler: "intact"|"split"|"uninformative"}``.
        min_split_votes: floor on the number of "split" votes required.

    Returns:
        dict with ``n_intact``, ``n_split``, ``n_uninformative``, ``n_informative``,
        ``verdict`` and ``reason``.  ``verdict`` is ``"mis_join_candidate"`` when a
        majority of informative assemblies split the contig, ``"corroborated"``
        when a majority reproduce it intact, and ``"unresolved"`` when there is
        too little information to say.

    Only informative votes count toward the majority, so an assembly that is
    too fragmented to have an opinion neither supports nor opposes the join.
    """
    n_intact = sum(1 for v in votes.values() if v == "intact")
    n_split = sum(1 for v in votes.values() if v == "split")
    n_uninf = sum(1 for v in votes.values() if v == "uninformative")
    n_inf = n_intact + n_split

    if n_inf == 0:
        verdict, reason = "unresolved", "no assembly produced an informative alignment"
    elif n_split >= min_split_votes and n_split > n_intact:
        verdict = "mis_join_candidate"
        reason = f"{n_split} of {n_inf} informative assemblies split this contig"
    elif n_intact >= n_split:
        verdict = "corroborated"
        reason = f"{n_intact} of {n_inf} informative assemblies reproduce it intact"
    else:
        verdict = "unresolved"
        reason = f"split by {n_split} of {n_inf} but below the {min_split_votes}-vote floor"

    return {"n_intact": n_intact, "n_split": n_split,
            "n_uninformative": n_uninf, "n_informative": n_inf,
            "verdict": verdict, "reason": reason}


def corroborated_t2t_count(t2t_contigs, verdicts, mode="exclude"):
    """Number of strict-T2T contigs that survive the concordance check.

    Args:
        t2t_contigs: iterable of contig names classified strict_t2t.
        verdicts: ``{contig: verdict_dict}`` from :func:`concordance_verdict`.
        mode: ``"exclude"`` discounts mis-join candidates from the count;
            ``"flag"`` reports them but leaves the count unchanged.

    Returns:
        ``(count, excluded)`` where *excluded* lists the discounted contigs.

    An ``unresolved`` contig is always counted, so a contig is only ever
    discounted on positive evidence of disagreement.
    """
    t2t = list(t2t_contigs)
    if mode != "exclude":
        return len(t2t), []
    excluded = [c for c in t2t
                if verdicts.get(c, {}).get("verdict") == "mis_join_candidate"]
    return len(t2t) - len(excluded), excluded


# ── Fusion detection against a --compare genome ──────────────────────────────

def detect_fusions(rows, min_query_frac=0.50, min_aligned_bp=100000):
    """Find contigs that absorb two or more sequences of the compare genome.

    ``contig_to_contig.tsv`` reports each compare sequence's best target, so a
    contig that receives two compare chromosomes is recorded as two independent
    "1-to-1" rows and the fusion is invisible.  This detects that pattern.

    Args:
        rows: iterable of dicts with ``compare_contig``, ``compare_len``,
            ``target_contig``, ``target_len`` and ``aligned_bp``.
        min_query_frac: a compare sequence must have this fraction of its own
            length aligned to the target to count as landing on it.
        min_aligned_bp: absolute floor, so short sequences and repeats do not
            create spurious fusion calls.

    Returns:
        list of dicts, one per suspect target contig, sorted by target name.
    """
    by_target = {}
    for r in rows:
        clen = int(r.get("compare_len") or 0)
        albp = int(r.get("aligned_bp") or 0)
        tgt = r.get("target_contig")
        if not tgt or clen <= 0:
            continue
        if albp < min_aligned_bp or albp < min_query_frac * clen:
            continue
        by_target.setdefault(tgt, []).append(r)

    out = []
    for tgt, members in sorted(by_target.items()):
        if len(members) < 2:
            continue
        tlen = int(members[0].get("target_len") or 0)
        summed = sum(int(m.get("compare_len") or 0) for m in members)
        out.append({
            "target_contig": tgt,
            "target_len": tlen,
            "n_compare_sequences": len(members),
            "compare_sequences": ",".join(str(m.get("compare_contig")) for m in members),
            "compare_total_len": summed,
            "length_excess_vs_sum": tlen - summed,
            # Junction positions implied by laying the compare sequences end to
            # end; the true breakpoint is near one of these offsets.
            "implied_junctions": ",".join(
                str(v) for v in _cumulative_offsets(
                    [int(m.get("compare_len") or 0) for m in members])),
            "reason": (f"{len(members)} compare sequences each align substantially "
                       f"to this single contig — possible mis-join"),
        })
    return out


def _cumulative_offsets(lengths):
    """Interior cumulative offsets, i.e. the implied junction coordinates."""
    offs = []
    run = 0
    for L in lengths[:-1]:
        run += L
        offs.append(run)
    return offs


# ── Contaminant / foreign-sequence screening ─────────────────────────────────

#: A contig whose median depth falls below this fraction of the assembly's
#: modal depth is anomalous — too shallow to be a chromosome of the sequenced
#: organism at the observed coverage.
DEPTH_RATIO_MAX = 0.34
#: Absolute GC deviation (percentage points) from the assembly mode that marks
#: a contig as compositionally foreign.
GC_DEV_MAX = 8.0
#: Contigs shorter than this are ignored, since short contigs are noisy in both
#: depth and GC.
MIN_SCREEN_LEN = 10000


def length_weighted_median(pairs):
    """Median of ``(value, weight)`` pairs, weighting by *weight*.

    Used for the assembly's modal depth and GC so that a few small aberrant
    contigs cannot move the baseline they are being compared against.
    """
    items = [(float(v), float(w)) for v, w in pairs if w and w > 0]
    if not items:
        return 0.0
    items.sort()
    half = sum(w for _, w in items) / 2.0
    run = 0.0
    for v, w in items:
        run += w
        if run >= half:
            return v
    return items[-1][0]


def screen_contaminants(contigs, depth_ratio_max=DEPTH_RATIO_MAX,
                        gc_dev_max=GC_DEV_MAX, min_len=MIN_SCREEN_LEN):
    """Flag contigs that look like foreign sequence rather than target genome.

    Args:
        contigs: iterable of dicts with ``contig``, ``length``, ``median_cov``,
            ``gc`` and optionally ``aligned_to_compare`` (bool) and
            ``busco_count`` (int).  Extra keys are ignored.
        depth_ratio_max / gc_dev_max / min_len: see module constants.

    Returns:
        ``(flagged, baseline)``.  *flagged* is a list of dicts with the
        triggering signals and a ``signals`` count; *baseline* records the modal
        depth and GC the comparison was made against.

    Read depth and GC are independent signals, so a contig is reported when
    either fires and the report states which.  Nothing is removed here: the
    caller decides, and TACO's convention is to write a filtered assembly
    alongside the unfiltered one rather than to delete sequence.
    """
    recs = [c for c in contigs if int(c.get("length") or 0) > 0]
    if not recs:
        return [], {"modal_depth": 0.0, "modal_gc": 0.0, "n_contigs": 0}

    # ``median_cov`` of None means "no coverage record for this contig", which
    # is not the same as zero coverage.  Treating a missing record as 0x would
    # flag it as a contaminant on the strength of absent data, so contigs
    # without depth are excluded from the baseline and only screened on GC.
    modal_depth = length_weighted_median(
        [(c.get("median_cov"), c.get("length") or 0) for c in recs
         if c.get("median_cov") is not None])
    modal_gc = length_weighted_median(
        [(c.get("gc") or 0, c.get("length") or 0) for c in recs])
    baseline = {"modal_depth": modal_depth, "modal_gc": modal_gc,
                "n_contigs": len(recs)}

    flagged = []
    for c in recs:
        length = int(c["length"])
        if length < min_len:
            continue
        raw_depth = c.get("median_cov")
        has_depth = raw_depth is not None
        depth = float(raw_depth) if has_depth else 0.0
        gc = float(c.get("gc") or 0)
        depth_informative = has_depth and modal_depth > 0
        ratio = (depth / modal_depth) if depth_informative else 1.0
        gc_dev = abs(gc - modal_gc)

        signals = []
        if depth_informative and ratio < depth_ratio_max:
            signals.append(f"depth {depth:.0f}x is {ratio:.2f} of modal {modal_depth:.0f}x")
        if gc_dev > gc_dev_max:
            signals.append(f"GC {gc:.1f}% deviates {gc_dev:.1f} points from modal {modal_gc:.1f}%")
        if c.get("aligned_to_compare") is False:
            signals.append("no alignment to the compare genome")
        if c.get("busco_count") is not None and int(c["busco_count"]) == 0 and length >= 1000000:
            signals.append(f"no BUSCO genes over {length:,} bp")

        # Depth or GC alone is enough to report; the other keys only corroborate.
        primary = (depth_informative and ratio < depth_ratio_max) or (gc_dev > gc_dev_max)
        if primary:
            flagged.append({
                "contig": c.get("contig"),
                "length": length,
                "median_cov": depth,
                "depth_ratio": round(ratio, 4),
                "gc": gc,
                "gc_deviation": round(gc_dev, 2),
                "n_signals": len(signals),
                "signals": "; ".join(signals),
            })
    flagged.sort(key=lambda r: -r["length"])
    return flagged, baseline
