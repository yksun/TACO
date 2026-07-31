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

Two independent checks are provided:

``concordance_verdict``
    Does another assembly split this contig?  Majority vote across assemblers.
``detect_fusions``
    In ``--compare`` mode, do two reference sequences land on one contig?

Contaminant screening and the machinery that *acts* on a mis-join — breakpoint
estimation and read-spanning confirmation — live in :mod:`taco.purify`.
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
#: Two "substantial" targets may overlap on the query by at most this fraction
#: of the smaller one before they are treated as a repeat rather than a split.
MAX_SEGMENT_OVERLAP_FRAC = 0.50


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
    if len(substantial) < 2:
        return "uninformative"

    # Two substantial targets are only evidence of a split if they occupy
    # *different* parts of the query.  Targets that pile up on the same region
    # are a repeat family, and calling that a split would break real contigs
    # wherever a segmental duplication happens to be present in two contigs of
    # another assembly.
    ranked = sorted(substantial, key=lambda t: -covers[t])[:2]
    a, b = target_intervals[ranked[0]], target_intervals[ranked[1]]
    overlap = _interval_overlap(a, b)
    if overlap >= MAX_SEGMENT_OVERLAP_FRAC * min(covers[ranked[0]],
                                                 covers[ranked[1]]):
        return "uninformative"
    return "split"


def _interval_overlap(ivs_a, ivs_b):
    """Bases covered by both interval sets."""
    if not ivs_a or not ivs_b:
        return 0
    total = 0
    for s1, e1 in ivs_a:
        for s2, e2 in ivs_b:
            lo, hi = max(s1, s2), min(e1, e2)
            if hi > lo:
                total += hi - lo
    return total


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

    A tie among informative voters counts as a flag, not as corroboration.  In
    v1.3.7 the majority test was strictly ``n_split > n_intact``, so an even
    split — two assemblies breaking the contig and two reproducing it — was
    silently reported as corroborated and the candidate escaped scrutiny
    entirely.  Flagging a tie is safe because being flagged only escalates a
    contig to investigation: :mod:`taco.purify` still requires an agreed
    breakpoint and read evidence before anything is broken.
    """
    n_intact = sum(1 for v in votes.values() if v == "intact")
    n_split = sum(1 for v in votes.values() if v == "split")
    n_uninf = sum(1 for v in votes.values() if v == "uninformative")
    n_inf = n_intact + n_split

    if n_inf == 0:
        verdict, reason = "unresolved", "no assembly produced an informative alignment"
    elif n_split >= min_split_votes and n_split >= n_intact:
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


# ── Contaminant screening moved to taco.purify in v1.4.0 ─────────────────────
#
# ``screen_contaminants`` and ``length_weighted_median`` lived here in v1.3.7.
# They now live in :mod:`taco.purify`, rewritten, because the versions here had
# two defects that only a rewrite could fix: a length-weighted median baseline
# inverts host and contaminant once foreign sequence passes half the assembly,
# and an OR over a single depth or GC signal at fixed absolute thresholds both
# condemns real atypical chromosomes and cannot adapt across genomes.
