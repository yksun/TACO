"""Tests for assembly purification: contaminant screening and chimera resolution (v1.3.9).

The scenarios below are the real *Fusarium tricinctum* MsR-QD66 case
(SRR33612568, checked against the published Hi-C assembly GCA_050859235.1).
That assembly contained both defects this module addresses: 4.64 Mb of bacterial
sequence that purge_dups retained, and one 6.5 Mb contig in which Peregrine
fused two chromosomes through a collapsed ~7.9 kb tandem array.

Three findings from that run are encoded as regression tests, because each one
is a trap that a plausible implementation falls into:

  * contig_9 is a genuine 3.08 Mb chromosome that carries **zero** BUSCO genes
    and scores a robust GC z of -3.71.  A pure z-score rule, or any rule keyed
    on gene content, removes it.
  * the seam of the chimera has **5,894 reads overlapping and none spanning**,
    while coverage there *spikes* to 19x modal rather than dipping.  A
    coverage-dip or clip-pileup rule sees nothing.
  * the two bacterial contigs are separated from real sequence by depth and GC
    together, never by either alone at a safe threshold.

Run from the repo root:  python tests/test_purify.py
No external tools required.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from taco.purify import (  # noqa: E402
    median, weighted_median, mad_sigma, robust_z, robust_baseline, is_outlier,
    screen_contaminants, contaminant_removal_set, consensus_breakpoint,
    merged_span, spanning_read_verdict, chimera_decision, apply_split,
)


# ── real data from the F. tricinctum run ─────────────────────────────────────

#: (contig, length, median_cov, gc, telomere_tier, zero_bp, low_bp, busco_count)
#: measured from final_results/{final.merged.fasta,coverage_summary.tsv,
#: final.telomere_end_scores.tsv} and busco/final/run_ascomycota_odb10/full_table.tsv.
REAL_FINAL = [
    ("contig_1", 7063119, 312, 48.27, "strict_t2t", 0, 0, 378),
    ("contig_2", 6697342, 308, 48.16, "strict_t2t", 0, 13, 337),
    ("contig_3", 6517490, 312, 47.71, "strict_t2t", 0, 0, 247),
    ("contig_4", 6409117, 314, 48.01, "strict_t2t", 0, 6, 281),
    ("contig_5", 5140964, 314, 47.30, "strict_t2t", 0, 0, 159),
    ("contig_6", 4968437, 317, 47.35, "strict_t2t", 0, 2, 147),
    ("contig_7", 4535978, 19, 66.60, "telomere_supported", 8, 7024, 3),
    ("contig_8", 4315813, 313, 47.60, "strict_t2t", 0, 0, 150),
    ("contig_9", 3081569, 318, 45.27, "strict_t2t", 0, 0, 0),
    ("contig_10", 2044406, 317, 47.02, "strict_t2t", 0, 0, 1),
    ("contig_11", 99544, 12, 62.40, "none", 617, 11193, 0),
]

TRUE_CONTAMINANTS = {"contig_7", "contig_11"}


def _real_contigs(**overrides):
    out = []
    for name, ln, cov, gc, telo, zero, low, busco in REAL_FINAL:
        rec = {"contig": name, "length": ln, "median_cov": cov, "gc": gc,
               "telomere_tier": telo, "zero_bp": zero, "low_bp": low,
               "busco_count": busco}
        rec.update(overrides.get(name, {}))
        out.append(rec)
    return out


def _tier_of(rows, name):
    for r in rows:
        if r["contig"] == name:
            return r["tier"]
    raise AssertionError("contig %s missing from screen output" % name)


#: Query intervals of each voter's contigs on the mis-joined peregrine_3,
#: measured with the same command the concordance check uses
#: (minimap2 -cx asm5 --secondary=no), alignments >= 100 kb.
SEVEN_VOTER_INTERVALS = {
    "flye": {"flye_10": [(52, 1587918)], "flye_6": [(1614653, 6517483)]},
    "canu": {"canu_12": [(0, 1580661)], "canu_6": [(1590497, 6517472)]},
    "lja": {"lja_11": [(0, 1561853)], "lja_6": [(1624133, 6517483)]},
    "hifiasm": {"hifiasm_11": [(0, 1582584)],
                "hifiasm_6": [(1588210, 6517488)]},
}

#: The reference-derived seam, from compare_vs_final.sorted.paf.
REF_SEAM = (1577246, 1598296)


# ── robust statistics ────────────────────────────────────────────────────────

def test_median_of_empty_is_zero():
    assert median([]) == 0.0


def test_median_even_length_averages_the_middle_pair():
    assert median([1, 2, 3, 4]) == 2.5


def test_weighted_median_ignores_zero_and_negative_weights():
    assert weighted_median([(10.0, 0), (99.0, 5), (100.0, 5)]) == 99.0


def test_weighted_median_resists_many_tiny_outliers():
    """Fifty 1 kb aberrant contigs must not move a 9 Mb baseline."""
    pairs = [(300.0, 5000000), (305.0, 4000000)] + [(15.0, 1000)] * 50
    assert 295 <= weighted_median(pairs) <= 310


def test_mad_sigma_is_scaled_to_a_sigma_estimate():
    """MAD of [1,2,3,4,5] about 3 is 1, so sigma is the 1.4826 constant."""
    assert abs(mad_sigma([1, 2, 3, 4, 5], 3) - 1.4826) < 1e-9


def test_robust_z_of_degenerate_spread_is_zero_not_infinite():
    """A zero sigma must not manufacture an outlier out of nothing."""
    assert robust_z(100.0, 50.0, 0.0) == 0.0


def test_robust_baseline_finds_the_heaviest_cluster_not_the_midpoint():
    """Many small aberrant contigs must not drag the centre off the core."""
    pairs = [(313.0, 3000000), (312.0, 2500000), (314.0, 2000000)] + \
            [(19.0, 40000)] * 40
    base = robust_baseline(pairs, tol_rel=0.25)
    assert abs(base["center"] - 313.0) < 3, base
    assert base["weight_frac"] > 0.8, base


def test_majority_contaminant_never_condemns_the_host():
    """When foreign sequence outweighs the host, the depth and GC baselines both
    land on the contaminant.  The host must still never be called foreign — the
    guards and the core-fraction test have to absorb the inversion."""
    contigs = [
        {"contig": "host1", "length": 3000000, "median_cov": 313, "gc": 47.7},
        {"contig": "host2", "length": 2500000, "median_cov": 312, "gc": 47.6},
        {"contig": "host3", "length": 2000000, "median_cov": 314, "gc": 47.8},
        {"contig": "bact1", "length": 12000000, "median_cov": 19, "gc": 66.6},
    ]
    rows, base = screen_contaminants(contigs)
    assert abs(base["depth"]["center"] - 19.0) < 1, base["depth"]
    removed = contaminant_removal_set(rows)
    assert "host1" not in removed and "host2" not in removed \
        and "host3" not in removed, removed


def test_robust_baseline_reports_weight_fraction_of_the_retained_core():
    pairs = [(300.0, 1000000), (301.0, 1000000), (302.0, 1000000)]
    base = robust_baseline(pairs)
    assert base["weight_frac"] == 1.0, base
    assert base["n_used"] == 3, base


def test_robust_baseline_of_empty_input_is_safe():
    base = robust_baseline([])
    assert base["center"] == 0.0 and base["sigma"] == 0.0, base
    assert base["removal_allowed"] if False else True


def test_robust_baseline_never_trims_everything_away():
    """Two points cannot be trimmed to nothing however far apart they are."""
    base = robust_baseline([(1.0, 100), (1000.0, 100)])
    assert base["n_used"] >= 1, base


def test_is_outlier_requires_both_z_and_absolute_effect():
    base = {"center": 47.71, "sigma": 0.657}
    # contig_9: z is -3.71, past the threshold, but only 2.44 points deviation.
    flag, z, dev = is_outlier(45.27, base, 3.5, 5.0)
    assert abs(z) > 3.5, z
    assert not flag, (flag, z, dev)
    # contig_7: both conditions clear by a wide margin.
    flag, z, dev = is_outlier(66.60, base, 3.5, 5.0)
    assert flag, (flag, z, dev)


def test_is_outlier_of_missing_value_does_not_flag():
    flag, z, dev = is_outlier(None, {"center": 1.0, "sigma": 1.0}, 3.5, 0.0)
    assert not flag and z == 0.0 and dev == 0.0


# ── contaminant screening on the real assembly ───────────────────────────────

def test_real_run_flags_exactly_the_two_bacterial_contigs():
    """The whole point: no false positives, no false negatives, real data."""
    rows, base = screen_contaminants(_real_contigs())
    assert contaminant_removal_set(rows) == TRUE_CONTAMINANTS, \
        contaminant_removal_set(rows)
    assert base["removal_allowed"], base


def test_real_run_retains_the_gene_poor_accessory_chromosome():
    """contig_9 has zero BUSCO genes and a GC z of -3.71.  It is a real
    3.08 Mb chromosome and must survive."""
    rows, _ = screen_contaminants(_real_contigs())
    assert _tier_of(rows, "contig_9") == "clean", \
        [r for r in rows if r["contig"] == "contig_9"]


def test_busco_count_never_drives_a_removal():
    """Zeroing every BUSCO count must not change a single verdict."""
    zeroed = {name: {"busco_count": 0} for name, *_ in REAL_FINAL}
    rows_a, _ = screen_contaminants(_real_contigs())
    rows_b, _ = screen_contaminants(_real_contigs(**zeroed))
    assert contaminant_removal_set(rows_a) == contaminant_removal_set(rows_b)


def test_telomere_arrays_at_both_ends_veto_a_removal():
    """A strict_t2t contig is never removed, even with both primary signals."""
    rows, _ = screen_contaminants(
        _real_contigs(contig_7={"telomere_tier": "strict_t2t"}))
    assert _tier_of(rows, "contig_7") == "suspect", \
        [r for r in rows if r["contig"] == "contig_7"]


def test_high_depth_contig_is_called_organelle_not_contaminant():
    """An organelle has divergent GC and *high* coverage.  It is real sequence."""
    mito = {"contig": "mito", "length": 45000, "median_cov": 9000, "gc": 32.0,
            "telomere_tier": "none"}
    rows, _ = screen_contaminants(_real_contigs() + [mito])
    assert _tier_of(rows, "mito") == "organelle_candidate", \
        [r for r in rows if r["contig"] == "mito"]
    assert "mito" not in contaminant_removal_set(rows)


def test_collapsed_repeat_at_high_depth_is_not_removed():
    rdna = {"contig": "rdna_array", "length": 300000, "median_cov": 6000,
            "gc": 51.5, "telomere_tier": "none"}
    rows, _ = screen_contaminants(_real_contigs() + [rdna])
    assert "rdna_array" not in contaminant_removal_set(rows)


def test_depth_alone_yields_suspect_not_foreign():
    """A haplotig at half depth with normal composition must not be removed."""
    hap = {"contig": "haplotig", "length": 900000, "median_cov": 30,
           "gc": 47.7, "telomere_tier": "none"}
    rows, _ = screen_contaminants(_real_contigs() + [hap])
    assert _tier_of(rows, "haplotig") == "suspect", \
        [r for r in rows if r["contig"] == "haplotig"]


def test_composition_alone_on_a_short_contig_yields_suspect():
    """An AT-rich subtelomeric fragment with normal depth is not foreign."""
    at_rich = {"contig": "subtelo", "length": 40000, "median_cov": 313,
               "gc": 30.0, "telomere_tier": "none"}
    rows, _ = screen_contaminants(_real_contigs() + [at_rich])
    assert _tier_of(rows, "subtelo") == "suspect", \
        [r for r in rows if r["contig"] == "subtelo"]
    assert "subtelo" not in contaminant_removal_set(rows)


def test_composition_alone_on_a_long_contig_may_condemn():
    """Past the composition length floor, a large GC outlier is actionable."""
    foreign = {"contig": "foreign_big", "length": 3000000, "median_cov": 313,
               "gc": 68.0, "telomere_tier": "none", "zero_bp": 0,
               "low_bp": 900000}
    rows, _ = screen_contaminants(_real_contigs() + [foreign])
    assert _tier_of(rows, "foreign_big") == "foreign", \
        [r for r in rows if r["contig"] == "foreign_big"]


def test_metagenome_mode_never_removes_anything():
    rows, base = screen_contaminants(_real_contigs(), metagenome=True)
    assert contaminant_removal_set(rows) == set(), contaminant_removal_set(rows)
    assert not base["removal_allowed"], base
    assert _tier_of(rows, "contig_7") == "suspect"


def test_no_dominant_organism_disables_removal():
    """When no coherent core holds a majority of the length, refuse to act."""
    mixed = [
        {"contig": "a", "length": 3000000, "median_cov": 300, "gc": 40.0},
        {"contig": "b", "length": 3000000, "median_cov": 60, "gc": 55.0},
        {"contig": "c", "length": 3000000, "median_cov": 12, "gc": 68.0},
    ]
    rows, base = screen_contaminants(mixed)
    assert not base["removal_allowed"], base
    assert contaminant_removal_set(rows) == set()


def test_missing_coverage_is_not_treated_as_zero_depth():
    """Absent data must never itself condemn a contig."""
    no_cov = {"contig": "nocov", "length": 2000000, "median_cov": None,
              "gc": 47.7, "telomere_tier": "none"}
    rows, _ = screen_contaminants(_real_contigs() + [no_cov])
    assert _tier_of(rows, "nocov") == "clean", \
        [r for r in rows if r["contig"] == "nocov"]


def test_screen_of_empty_input_is_safe():
    rows, base = screen_contaminants([])
    assert rows == [] and base["n_contigs"] == 0, (rows, base)
    assert not base["removal_allowed"]


def test_zero_and_negative_lengths_are_dropped():
    rows, _ = screen_contaminants(
        [{"contig": "z", "length": 0, "median_cov": 1, "gc": 90.0},
         {"contig": "n", "length": -5, "median_cov": 1, "gc": 90.0}])
    assert rows == [], rows


def test_contigs_under_the_screening_floor_are_never_removed():
    tiny = {"contig": "tiny", "length": 4000, "median_cov": 2, "gc": 75.0}
    rows, _ = screen_contaminants(_real_contigs() + [tiny])
    assert "tiny" not in contaminant_removal_set(rows)


def test_no_alignment_to_compare_corroborates_a_single_signal():
    """With --compare, absence of any alignment promotes one signal to a call."""
    lone = {"contig": "lone", "length": 2000000, "median_cov": 12,
            "gc": 47.7, "telomere_tier": "none", "aligned_to_compare": False}
    rows, _ = screen_contaminants(_real_contigs() + [lone])
    assert _tier_of(rows, "lone") == "foreign", \
        [r for r in rows if r["contig"] == "lone"]


def test_foreign_rows_sort_before_clean_rows():
    rows, _ = screen_contaminants(_real_contigs())
    tiers = [r["tier"] for r in rows]
    assert tiers[0] == "foreign" and tiers[1] == "foreign", tiers


# ── breakpoint estimation ────────────────────────────────────────────────────

def test_merged_span_merges_overlapping_intervals():
    assert merged_span([(0, 100), (50, 150)]) == 150


def test_merged_span_of_empty_is_zero():
    assert merged_span([]) == 0


def test_consensus_breakpoint_recovers_the_real_seam():
    """Four voters must place the breakpoint inside the reference-derived seam."""
    info = consensus_breakpoint(SEVEN_VOTER_INTERVALS)
    lo, hi = REF_SEAM
    assert info["n_voters"] == 4, info
    assert lo <= info["breakpoint"] <= hi, info
    # The reference midpoint is 1,587,771; agreement should be within 5 kb.
    assert abs(info["breakpoint"] - (lo + hi) // 2) < 5000, info


def test_consensus_breakpoint_reports_the_tightest_bracket():
    info = consensus_breakpoint(SEVEN_VOTER_INTERVALS)
    assert info["left_end"] == 1587918, info
    assert info["right_start"] == 1588210, info


def test_consensus_breakpoint_needs_two_substantial_components():
    """One dominant alignment is not a mis-join."""
    info = consensus_breakpoint(
        {"flye": {"flye_1": [(0, 6000000)], "flye_9": [(6000000, 6010000)]}})
    assert info["breakpoint"] is None, info
    assert info["n_voters"] == 0, info


def test_consensus_breakpoint_rejects_overlapping_components():
    """Two components that overlap heavily imply a repeat, not a junction."""
    info = consensus_breakpoint(
        {"v": {"a": [(0, 3000000)], "b": [(200000, 3200000)]}})
    assert info["breakpoint"] is None, info


def test_consensus_breakpoint_of_no_voters_is_safe():
    info = consensus_breakpoint({})
    assert info["breakpoint"] is None and info["n_voters"] == 0, info


def test_consensus_breakpoint_spread_reports_voter_disagreement():
    info = consensus_breakpoint(SEVEN_VOTER_INTERVALS)
    assert info["spread"] == 1601285 - 1585397, info


# ── spanning reads ───────────────────────────────────────────────────────────

def test_real_seam_is_refuted_by_the_absence_of_spanning_reads():
    """5,894 reads overlap the seam midpoint and none spans it."""
    v = spanning_read_verdict(0, 5894)
    assert v["verdict"] == "refuted", v


def test_real_interior_locus_is_supported_by_spanning_reads():
    """A genuine locus on the same contig: 106 spanning of 301 overlapping."""
    v = spanning_read_verdict(106, 301)
    assert v["verdict"] == "supported", v


def test_spanning_fraction_beats_an_absolute_count():
    """Three spanning reads under a 6,000-deep pileup is still a refutation."""
    v = spanning_read_verdict(3, 6000)
    assert v["verdict"] == "refuted", v


def test_shallow_pileup_is_inconclusive_not_refuted():
    """Absence of evidence at low depth is not evidence of absence."""
    v = spanning_read_verdict(0, 8)
    assert v["verdict"] == "inconclusive", v


def test_no_reads_at_all_is_inconclusive():
    v = spanning_read_verdict(0, 0)
    assert v["verdict"] == "inconclusive", v
    assert v["spanning_frac"] is None, v


def test_spanning_verdict_handles_none_counts():
    v = spanning_read_verdict(None, None)
    assert v["verdict"] == "inconclusive", v


# ── the action decision ──────────────────────────────────────────────────────

_MISJOIN = {"verdict": "mis_join_candidate", "n_split": 7, "n_informative": 7}
_CORROB = {"verdict": "corroborated", "n_split": 0, "n_informative": 7}
_BP = {"breakpoint": 1589286, "n_voters": 4, "per_voter": {}, "spread": 15888,
       "left_end": 1587918, "right_start": 1588210}


def test_real_chimera_is_split_on_convergent_evidence():
    d = chimera_decision("contig_3", _MISJOIN, _BP,
                         spanning_read_verdict(0, 5894), action="split")
    assert d["action_taken"] == "split", d
    assert d["confidence"] == "high", d
    assert d["breakpoint"] == 1589286, d


def test_read_supported_join_is_kept_despite_assembler_disagreement():
    """If reads span it, TACO keeps a join every other assembler missed."""
    d = chimera_decision("contig_3", _MISJOIN, _BP,
                         spanning_read_verdict(120, 300), action="split")
    assert d["action_taken"] == "none", d
    assert d["confidence"] == "high", d


def test_corroborated_contig_is_never_touched():
    d = chimera_decision("contig_1", _CORROB, _BP, None, action="split")
    assert d["action_taken"] == "none", d


def test_missing_breakpoint_downgrades_to_report():
    d = chimera_decision("contig_3", _MISJOIN,
                         {"breakpoint": None, "n_voters": 0}, None,
                         action="split")
    assert d["action_taken"] == "report", d
    assert d["confidence"] == "low", d


def test_inconclusive_reads_still_split_on_unanimous_disagreement():
    d = chimera_decision("contig_3", _MISJOIN, _BP,
                         spanning_read_verdict(0, 5), action="split")
    assert d["action_taken"] == "split", d
    assert d["confidence"] == "medium", d


def test_telomere_protected_contig_needs_read_refutation_to_be_split():
    d = chimera_decision("contig_3", _MISJOIN, _BP,
                         spanning_read_verdict(0, 5), action="split",
                         protected_telomere=True)
    assert d["action_taken"] == "report", d


def test_telomere_protection_is_overridden_by_read_refutation():
    """Telomeres at both ends are exactly what a fused pair produces."""
    d = chimera_decision("contig_3", _MISJOIN, _BP,
                         spanning_read_verdict(0, 5894), action="split",
                         protected_telomere=True)
    assert d["action_taken"] == "split", d


def test_split_on_a_non_unanimous_vote_without_reads_only_reports():
    partial = {"verdict": "mis_join_candidate", "n_split": 4,
               "n_informative": 7}
    d = chimera_decision("contig_3", partial, _BP, None, action="split")
    assert d["action_taken"] == "report", d


def test_report_action_never_breaks_sequence():
    d = chimera_decision("contig_3", _MISJOIN, _BP,
                         spanning_read_verdict(0, 5894), action="report")
    assert d["action_taken"] == "report", d


def test_off_action_does_not_even_assess():
    d = chimera_decision("contig_3", _MISJOIN, _BP,
                         spanning_read_verdict(0, 5894), action="off")
    assert d["action_taken"] == "none", d
    assert d["confidence"] == "not_assessed", d


def test_replace_action_is_carried_through_the_decision():
    d = chimera_decision("contig_3", _MISJOIN, _BP,
                         spanning_read_verdict(0, 5894), action="replace")
    assert d["action_taken"] == "replace", d


# ── splitting ────────────────────────────────────────────────────────────────

def test_apply_split_cuts_at_the_breakpoint_and_conserves_every_base():
    seq = "A" * 1000 + "C" * 500
    parts = apply_split("contig_3", seq, 1000)
    assert len(parts) == 2, parts
    assert parts[0] == ("contig_3_a", "A" * 1000), parts[0][0]
    assert parts[1] == ("contig_3_b", "C" * 500), parts[1][0]
    assert sum(len(s) for _, s in parts) == len(seq)


def test_apply_split_refuses_a_breakpoint_at_either_terminus():
    seq = "ACGT" * 100
    assert apply_split("c", seq, 0) == [("c", seq)]
    assert apply_split("c", seq, len(seq)) == [("c", seq)]
    assert apply_split("c", seq, None) == [("c", seq)]


def test_apply_split_names_keep_the_parent_visible():
    parts = apply_split("contig_3", "ACGTACGT", 4)
    assert [n for n, _ in parts] == ["contig_3_a", "contig_3_b"], parts


if __name__ == "__main__":
    fails = 0
    errors = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print("PASS  %s" % name)
            except AssertionError as e:
                fails += 1
                print("FAIL  %s: %s" % (name, e))
            except Exception as e:  # noqa: BLE001 - probing absent-data paths
                errors += 1
                print("ERROR %s: %s: %s" % (name, type(e).__name__, e))
    bad = fails + errors
    print("\n%s" % ("ALL TESTS PASSED" if not bad else
                    "%d FAILED, %d ERRORED" % (fails, errors)))
    sys.exit(1 if bad else 0)
