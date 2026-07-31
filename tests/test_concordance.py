"""Tests for cross-assembler concordance, mis-join and contaminant screening (v1.3.7).

The scenarios below are the real *Fusarium tricinctum* MsR-QD66 case
(SRR33612568, compared against GCA_050859235.1), in which the Peregrine
assembly presented one contig as telomere-to-telomere while every other
assembler split it in two.  Keeping these as regression tests means the
detection cannot silently break.

Run from the repo root:  python tests/test_concordance.py
No external tools required.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from taco.concordance import (  # noqa: E402
    merged_span, split_vote, concordance_verdict, corroborated_t2t_count,
    detect_fusions,
)


# ── interval arithmetic ──────────────────────────────────────────────────────

def test_merged_span_merges_overlaps():
    assert merged_span([]) == 0
    assert merged_span([(0, 100)]) == 100
    assert merged_span([(0, 100), (50, 150)]) == 150      # overlapping
    assert merged_span([(0, 100), (200, 300)]) == 200     # disjoint
    assert merged_span([(0, 100), (100, 200)]) == 200     # abutting


# ── per-assembly split voting ────────────────────────────────────────────────

def test_split_vote_detects_a_two_piece_split():
    """Peregrine contig_3 (6,517,490 bp) is two contigs in canu."""
    iv = {"canu_6": [(0, 4919096)], "canu_11": [(4930000, 6507246)]}
    assert split_vote(6517490, iv) == "split"


def test_split_vote_accepts_an_intact_chromosome():
    assert split_vote(7063119, {"canu_1": [(0, 7063110)]}) == "intact"


def test_split_vote_ignores_scattered_repeats():
    """Several small hits are not a structural split — they are repeats."""
    iv = {"x_1": [(0, 90000)], "x_2": [(100000, 180000)], "x_3": [(200000, 260000)]}
    assert split_vote(6000000, iv) == "uninformative"


def test_split_vote_needs_enough_coverage_to_have_an_opinion():
    # only 30% of the query aligns anywhere: no opinion
    assert split_vote(1000000, {"a": [(0, 300000)]}) == "uninformative"


def test_split_vote_handles_degenerate_input():
    assert split_vote(0, {"a": [(0, 10)]}) == "uninformative"
    assert split_vote(1000, {}) == "uninformative"


# ── combining votes ──────────────────────────────────────────────────────────

def test_unanimous_split_is_a_mis_join_candidate():
    """The real result: 7 assemblers split peregrine_3, none reproduced it."""
    votes = {a: "split" for a in
             ["canu", "flye", "hifiasm", "ipa", "lja", "nextDenovo", "raven"]}
    v = concordance_verdict(votes)
    assert v["verdict"] == "mis_join_candidate"
    assert v["n_split"] == 7 and v["n_intact"] == 0


def test_unanimous_intact_is_corroborated():
    v = concordance_verdict({a: "intact" for a in ["canu", "flye", "lja"]})
    assert v["verdict"] == "corroborated"


def test_single_dissenting_vote_does_not_flag():
    """peregrine_9: 3 intact, 1 split -> must stay corroborated.

    A lone disagreement is far more likely to be a fragmented assembly than a
    chimera, so the two-vote floor and the majority rule must both hold.
    """
    v = concordance_verdict({"canu": "intact", "flye": "intact",
                             "lja": "intact", "nextDenovo": "split"})
    assert v["verdict"] == "corroborated", v


def test_uninformative_votes_do_not_count_toward_the_majority():
    v = concordance_verdict({"a": "split", "b": "split",
                             "c": "uninformative", "d": "uninformative"})
    assert v["verdict"] == "mis_join_candidate"
    assert v["n_informative"] == 2


def test_no_informative_votes_is_unresolved_not_a_flag():
    """Fail safe: if alignment produced nothing, never penalize a contig."""
    v = concordance_verdict({a: "uninformative" for a in ["x", "y", "z"]})
    assert v["verdict"] == "unresolved"


def test_two_vote_floor_blocks_a_lone_split():
    v = concordance_verdict({"a": "split"})
    assert v["verdict"] != "mis_join_candidate"


# ── corrected T2T counts ─────────────────────────────────────────────────────

def test_corroborated_count_discounts_only_flagged_contigs():
    t2t = ["peregrine_%d" % i for i in (1, 2, 3, 4, 5, 6, 8, 9, 10)]
    verdicts = {"peregrine_3": {"verdict": "mis_join_candidate"}}
    n, excluded = corroborated_t2t_count(t2t, verdicts, mode="exclude")
    assert n == 8 and excluded == ["peregrine_3"]


def test_flag_mode_leaves_the_count_untouched():
    t2t = ["c1", "c2", "c3"]
    verdicts = {"c2": {"verdict": "mis_join_candidate"}}
    assert corroborated_t2t_count(t2t, verdicts, mode="flag") == (3, [])


def test_unresolved_contigs_are_still_counted():
    t2t = ["c1", "c2"]
    verdicts = {"c2": {"verdict": "unresolved"}}
    n, excluded = corroborated_t2t_count(t2t, verdicts, mode="exclude")
    assert n == 2 and excluded == []


# ── fusion detection against a compare genome ────────────────────────────────

def test_detect_fusions_finds_two_chromosomes_on_one_contig():
    """Reference chr6 and chr10 both land on contig_3."""
    rows = [
        {"compare_contig": "CM116343.1", "compare_len": 4919202,
         "target_contig": "contig_3", "target_len": 6517490, "aligned_bp": 4919096},
        {"compare_contig": "CM116347.1", "compare_len": 1577254,
         "target_contig": "contig_3", "target_len": 6517490, "aligned_bp": 1577246},
        {"compare_contig": "CM116338.1", "compare_len": 7063146,
         "target_contig": "contig_1", "target_len": 7063119, "aligned_bp": 7063120},
    ]
    out = detect_fusions(rows)
    assert len(out) == 1
    f = out[0]
    assert f["target_contig"] == "contig_3"
    assert f["n_compare_sequences"] == 2
    assert f["length_excess_vs_sum"] == 21034
    assert f["implied_junctions"] == "4919202"


def test_detect_fusions_ignores_partial_and_short_hits():
    """A sequence that only partly aligns must not create a fusion call."""
    rows = [
        {"compare_contig": "A", "compare_len": 5000000, "target_contig": "t1",
         "target_len": 5000000, "aligned_bp": 4900000},
        {"compare_contig": "B", "compare_len": 4000000, "target_contig": "t1",
         "target_len": 5000000, "aligned_bp": 150000},   # only 3.8% of B
    ]
    assert detect_fusions(rows) == []


def test_detect_fusions_returns_nothing_for_clean_one_to_one():
    rows = [{"compare_contig": "A", "compare_len": 1000000, "target_contig": "t1",
             "target_len": 1000000, "aligned_bp": 999000}]
    assert detect_fusions(rows) == []


# ── contaminant screening ────────────────────────────────────────────────────

REAL_FINAL = [
    # (contig, length, median_cov, gc) — the real 11 final contigs
    ("contig_1", 7063119, 312, 48.27), ("contig_2", 6697342, 308, 48.16),
    ("contig_3", 6517490, 312, 47.71), ("contig_4", 6409117, 314, 48.01),
    ("contig_5", 5140964, 314, 47.30), ("contig_6", 4968437, 317, 47.35),
    ("contig_7", 4535978, 19, 66.60), ("contig_8", 4315813, 313, 47.60),
    ("contig_9", 3081569, 318, 45.27), ("contig_10", 2044406, 317, 47.02),
    ("contig_11", 99544, 12, 62.40),
]




if __name__ == "__main__":
    fails = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print("PASS  %s" % name)
            except AssertionError as e:
                fails += 1
                print("FAIL  %s: %s" % (name, e))
    print("\n%s" % ("ALL TESTS PASSED" if not fails else "%d TEST(S) FAILED" % fails))
    sys.exit(1 if fails else 0)
