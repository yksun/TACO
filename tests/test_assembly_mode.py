"""Tests for --assembly-mode {both,primary,full} (v1.5.0).

TACO delivers two different genome representations. ``primary`` is a
nonredundant representative genome; ``full`` retains divergent
alternate/haplotype sequence rather than collapsing it. ``both`` is the default
and emits each of them from one set of assemblies.

Neither mode phases anything. Nothing here should be read as a haplotype
assignment.

Three properties carry the design and are asserted directly:

1. Selection scores what refinement CANNOT change above what it demonstrably
   fixes. Measured on the real run, backbone -> delivered: BUSCO C and Merqury
   completeness only degraded, while QV, T2T and duplication all improved.
2. Merqury completeness means opposite things in the two modes -- the objective
   itself in ``full``, evidence of retained redundancy in ``primary`` -- so it is
   weighted accordingly rather than treated as a neutral quality metric.
3. ``full`` never removes sequence for merely resembling a second haplotype,
   while contaminant screening and chimera resolution stay on in every mode.

Run from the repo root:  python tests/test_assembly_mode.py
No external tools required.
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from taco import policy  # noqa: E402
from taco.backbone import select_backbone, _compute_score  # noqa: E402
from taco.cli import parse_args  # noqa: E402
from taco.policy import (  # noqa: E402
    AssemblyPolicy, DEFAULT_MODE, DELIVERABLE_MODES, LEGACY_PRIMARY_WEIGHTS,
    MODE_BOTH, SCORING_PROFILES, VALID_MODES, haplotype_retention_advisory,
    resolve_weights, score_assembly, size_expectation, telomere_reward_cap,
)

G = 127_000_000
#: The two published NCBI assemblies used as ground truth below.
REF_HAPLOID = 122_318_233     # GCA_026914185.2, NCBI reference genome, 18 chrs
REF_DIKARYON = 253_538_837    # GCA_019358815.1, Pt76, 36 chrs = 2 x 18


# ── CLI ──────────────────────────────────────────────────────────────────────

BASE_ARGV = ["-g", "127m", "-t", "8", "--fastq", "reads.fq"]


def _parse(extra=()):
    saved = sys.argv
    try:
        sys.argv = ["TACO"] + BASE_ARGV + list(extra)
        return parse_args()
    finally:
        sys.argv = saved


def test_cli_defaults_to_both():
    """Neither objective is knowable from the reads, so TACO delivers both."""
    assert _parse().assembly_mode == MODE_BOTH
    assert DEFAULT_MODE == MODE_BOTH


def test_cli_accepts_each_mode_explicitly():
    for mode in ("both", "primary", "full"):
        assert _parse(["--assembly-mode", mode]).assembly_mode == mode


def test_cli_rejects_an_unknown_mode():
    for bad in ("phased", "hap1", "smart"):
        try:
            _parse(["--assembly-mode", bad])
        except SystemExit as e:
            assert e.code == 2, f"expected argparse exit 2 for {bad}, got {e.code}"
        else:
            raise AssertionError(f"--assembly-mode {bad} was accepted")


def test_no_purge_dups_still_exists():
    assert _parse(["--no-purge-dups"]).no_purge_dups is True


# ── 'both' is a driver, not a scoring profile ────────────────────────────────

def test_both_names_no_weight_set():
    """The two objectives disagree, so 'both' cannot resolve to one formula."""
    assert MODE_BOTH not in SCORING_PROFILES
    assert set(SCORING_PROFILES) == set(DELIVERABLE_MODES)
    try:
        resolve_weights(MODE_BOTH)
    except ValueError as e:
        assert "no single weight set" in str(e)
    else:
        raise AssertionError("resolve_weights('both') should refuse")


def test_valid_modes_are_the_deliverables_plus_both():
    assert set(VALID_MODES) == set(DELIVERABLE_MODES) | {MODE_BOTH}


def test_both_expands_to_one_policy_per_deliverable():
    subs = AssemblyPolicy(mode=MODE_BOTH, taxon="fungal").sub_policies()
    assert [p.mode for p in subs] == list(DELIVERABLE_MODES)
    assert all(not p.is_both for p in subs)


def test_a_single_mode_policy_is_its_own_only_deliverable():
    for mode in DELIVERABLE_MODES:
        assert [p.mode for p in AssemblyPolicy(mode=mode).sub_policies()] == [mode]


def test_both_refuses_to_answer_removal_questions():
    """The modes disagree about what may be removed, so 'both' must not guess."""
    p = AssemblyPolicy(mode=MODE_BOTH)
    for attr in ("purge_dups_enabled", "collapse_redundancy_enabled"):
        try:
            getattr(p, attr)
        except ValueError as e:
            assert "not defined for mode 'both'" in str(e)
        else:
            raise AssertionError(f"'both' answered {attr}")


def test_both_still_screens_contaminants_and_chimeras():
    p = AssemblyPolicy(mode=MODE_BOTH)
    assert p.contamination_screen_enabled
    assert p.chimera_resolution_enabled


def test_both_describes_every_deliverable():
    text = AssemblyPolicy(mode=MODE_BOTH).describe()
    assert text.startswith(f"Assembly representation mode: {MODE_BOTH}")
    for mode in DELIVERABLE_MODES:
        assert mode in text


def test_pipeline_splits_shared_from_per_mode_steps():
    """Steps 0-10 describe the reads; only 11-14 depend on the objective."""
    from taco.pipeline import PipelineRunner as P
    assert P.SHARED_STEPS == frozenset(range(0, 11))
    assert P.PER_MODE_STEPS == frozenset(range(11, 15))
    assert not (P.SHARED_STEPS & P.PER_MODE_STEPS)
    assert (P.SHARED_STEPS | P.PER_MODE_STEPS) == frozenset(range(0, 15))


# ── refinement-aware selection ───────────────────────────────────────────────
#
# Measured on the real Puccinia triticina run (backbone ipa -> delivered):
#   UNFIXABLE, only degraded : BUSCO C 96.3->96.2, completeness 91.74->91.51
#   FIXABLE, improved        : QV 40.50->47.85, T2T 1->4, BUSCO D 6.9->5.3

def test_primary_rewards_busco_complete_not_single():
    """BUSCO S is partly a function of duplication, which purge_dups removes,
    so rewarding S pre-purge rewards work already done."""
    for mode in DELIVERABLE_MODES:
        assert resolve_weights(mode, "fungal")["busco_term"] == policy.BUSCO_TERM_COMPLETE


def test_contiguity_outweighs_duplication_in_primary():
    """Refinement removes duplication; it cannot create contiguity."""
    w = resolve_weights("primary", "fungal")
    assert w["w_n50"] > w["w_busco_d"], (w["w_n50"], w["w_busco_d"])
    assert w["w_contigs"] < w["w_n50"]


def test_duplication_penalty_is_far_below_the_legacy_level():
    """It was 600/point -- the heaviest term after BUSCO -- for a defect that
    step 12H removes anyway."""
    assert LEGACY_PRIMARY_WEIGHTS["w_busco_d"] == 500
    assert resolve_weights("primary", "fungal")["w_busco_d"] <= 200


def test_qv_weight_reflects_that_it_is_a_log_scale():
    """A 31-point QV spread is a ~1300x difference in consensus error, which
    the legacy weight of 20/point effectively ignored."""
    assert LEGACY_PRIMARY_WEIGHTS["w_merqury_qv"] == 20
    for mode in DELIVERABLE_MODES:
        assert resolve_weights(mode, "fungal")["w_merqury_qv"] >= 100


def test_telomere_reward_is_discounted_because_refinement_supplies_it():
    """Telomere-pool rescue took the delivered assembly from 1 T2T to 4."""
    for mode in DELIVERABLE_MODES:
        w = resolve_weights(mode, "fungal")
        assert w["w_t2t"] < LEGACY_PRIMARY_WEIGHTS["w_t2t"]
        assert w["w_single_tel"] < LEGACY_PRIMARY_WEIGHTS["w_single_tel"]


def test_completeness_means_opposite_things_in_the_two_modes():
    """It is the fraction of READ k-mers present. A collapsed haploid
    representation of a heterozygous sample cannot hold both haplotypes', so a
    high value is the objective in full and redundancy in primary."""
    p = resolve_weights("primary", "fungal")["w_merqury_comp"]
    f = resolve_weights("full", "fungal")["w_merqury_comp"]
    assert f > 10 * p, (p, f)


def test_qv_is_kept_separate_from_completeness():
    for mode in DELIVERABLE_MODES:
        w = resolve_weights(mode, "fungal")
        assert w["w_merqury_comp"] != w["w_merqury_qv"]


def test_primary_selection_is_deliberately_not_v141_equivalent():
    """v1.5.0 changes how a backbone is chosen. The legacy weights are kept as a
    record of what the published v1.4.x runs used, and must not be resurrected
    silently by any profile."""
    w = resolve_weights("primary", "fungal")
    differing = [k for k in ("w_busco_d", "w_merqury_comp", "w_merqury_qv",
                             "w_n50", "w_t2t", "w_single_tel", "busco_term")
                 if w[k] != LEGACY_PRIMARY_WEIGHTS[k]]
    assert len(differing) >= 6, differing
    assert LEGACY_PRIMARY_WEIGHTS["busco_term"] == policy.BUSCO_TERM_SINGLE


def test_taxon_override_cannot_reimpose_a_duplication_penalty_on_full():
    assert resolve_weights("full", "fungal")["w_busco_d"] == 0
    assert resolve_weights("primary", "fungal")["w_busco_d"] > 0


def test_taxon_override_cannot_reimpose_an_absolute_contig_penalty_on_full():
    w = resolve_weights("full", "plant")
    assert w["w_contigs"] == 0 and w["w_contig_density"] > 0
    assert resolve_weights("primary", "plant")["w_contigs"] > 0


def test_taxon_presets_keep_their_relative_ordering():
    """Fungal strictest on duplication, plant most permissive (polyploidy)."""
    d = {t: resolve_weights("primary", t)["w_busco_d"]
         for t in ("fungal", "vertebrate", "plant")}
    assert d["fungal"] > d["vertebrate"] > d["plant"], d


def test_unknown_mode_raises():
    for bad in ("phased", "hap1", "", None):
        try:
            resolve_weights(bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"resolve_weights accepted {bad!r}")


# ── the real Puccinia triticina candidates ───────────────────────────────────
#
# Nine selectable assemblies from a completed run (-g 127m, --taxon fungal,
# --busco basidiomycota_odb10). The published references are above.

PT = {
 "ipa":        dict(busco_c=96.3, busco_s=89.5, busco_d=6.9,  contigs=329,   n50=2280396, t2t=1,  single_tel=48,   merqury_qv=40.4993, merqury_comp=91.7373, total_len=133691058),
 "nextDenovo": dict(busco_c=96.4, busco_s=87.0, busco_d=9.4,  contigs=241,   n50=1919314, t2t=0,  single_tel=27,   merqury_qv=56.5211, merqury_comp=93.8407, total_len=134027129),
 "hifiasm":    dict(busco_c=96.6, busco_s=91.4, busco_d=5.2,  contigs=519,   n50=7217872, t2t=10, single_tel=15,   merqury_qv=54.3543, merqury_comp=92.1773, total_len=147242165),
 "peregrine":  dict(busco_c=95.6, busco_s=78.8, busco_d=16.8, contigs=99,    n50=4075283, t2t=4,  single_tel=28,   merqury_qv=71.7249, merqury_comp=93.1071, total_len=147196254),
 "flye":       dict(busco_c=96.3, busco_s=31.7, busco_d=64.5, contigs=1068,  n50=495893,  t2t=0,  single_tel=66,   merqury_qv=56.9858, merqury_comp=98.7189, total_len=215910115),
 "lja":        dict(busco_c=96.7, busco_s=0.8,  busco_d=95.9, contigs=972,   n50=2609728, t2t=3,  single_tel=153,  merqury_qv=47.4854, merqury_comp=99.1194, total_len=253991215),
 "canu":       dict(busco_c=96.7, busco_s=5.5,  busco_d=91.2, contigs=1574,  n50=2547273, t2t=7,  single_tel=92,   merqury_qv=67.9531, merqury_comp=99.1841, total_len=271708550),
 "raven":      dict(busco_c=85.8, busco_s=82.1, busco_d=3.7,  contigs=2198,  n50=46413,   t2t=0,  single_tel=33,   merqury_qv=26.3035, merqury_comp=68.0004, total_len=86388243),
 "mbg":        dict(busco_c=85.3, busco_s=26.0, busco_d=59.2, contigs=46812, n50=13298,   t2t=26, single_tel=1668, merqury_qv=51.5389, merqury_comp=99.2552, total_len=287169572),
}
#: Candidates that collapsed the haplotypes. A full-mode win by one of these
#: would defeat the mode.
COLLAPSED = {"ipa", "nextDenovo", "hifiasm", "peregrine", "raven"}
#: Candidates that retained alternate sequence (BUSCO D >= 60, completeness >= 98).
RETAINING = {"lja", "canu", "flye"}


def _rank(mode):
    sc = {k: score_assembly(m, mode=mode, taxon="fungal", expected_haploid=G)[0]
          for k, m in PT.items()}
    return sorted(sc, key=sc.get, reverse=True), sc


def test_primary_selects_the_contiguous_accurate_assembly():
    """peregrine: 99 contigs, N50 4.08 Mb, QV 71.7. Its duplication (16.8%) is
    what purge_dups exists to remove."""
    order, _ = _rank("primary")
    assert order[0] == "peregrine", order


def test_primary_no_longer_selects_on_pre_purge_duplication():
    """The legacy weights chose ipa: 329 contigs, N50 2.28 Mb, QV 40.5 --
    delivering an assembly less contiguous and 24 QV points worse than the
    candidate they passed over."""
    order, _ = _rank("primary")
    assert order.index("peregrine") < order.index("ipa")
    assert PT["peregrine"]["contigs"] < PT["ipa"]["contigs"]
    assert PT["peregrine"]["merqury_qv"] > PT["ipa"]["merqury_qv"]


def test_primary_keeps_haplotype_retaining_assemblies_out_of_the_top():
    order, _ = _rank("primary")
    assert not (set(order[:4]) & RETAINING), order[:4]


def test_full_selects_the_haplotype_retaining_assembly():
    order, _ = _rank("full")
    assert order[0] in RETAINING, order
    assert order[0] == "lja", order


def test_full_matches_the_published_dikaryon_reference_size():
    order, _ = _rank("full")
    ratio = PT[order[0]]["total_len"] / REF_DIKARYON
    assert 0.97 <= ratio <= 1.03, ratio


def test_full_does_not_select_a_collapsed_assembly():
    """Without a dominant completeness term full mode has no positive signal for
    retention at all -- BUSCO C saturates near 96% and duplication and size are
    deliberately unscored -- and it picks a collapsed 147 Mb assembly."""
    order, sc = _rank("full")
    assert order[0] not in COLLAPSED, order
    best_collapsed = next(k for k in order if k in COLLAPSED)
    margin = sc[order[0]] - sc[best_collapsed]
    assert margin > 0.01 * abs(sc[order[0]]), (
        f"margin over {best_collapsed} is only {margin:,.0f}; too fragile")


def test_the_two_modes_disagree_on_the_same_candidates():
    assert _rank("primary")[0][0] != _rank("full")[0][0]


def test_the_fragmented_artifact_loses_in_every_mode():
    """It has the HIGHEST Merqury completeness of all nine (99.2552), so only
    fragmentation can separate it."""
    for mode in DELIVERABLE_MODES:
        order, sc = _rank(mode)
        assert order[-1] == "mbg", (mode, order)
    assert PT["mbg"]["merqury_comp"] > PT["lja"]["merqury_comp"]


def test_full_mode_treats_duplication_as_neutral():
    base = dict(busco_c=96.0, busco_s=96.0, busco_d=0.0, merqury_comp=95.0,
                merqury_qv=50.0, contigs=300, n50=2_500_000,
                total_len=130_000_000)
    dup = dict(base, busco_s=48.0, busco_d=48.0)    # same C, half duplicated
    a, _ = score_assembly(base, mode="full", taxon="fungal", expected_haploid=G)
    b, _ = score_assembly(dup, mode="full", taxon="fungal", expected_haploid=G)
    assert abs(a - b) < 1e-9


def test_busco_c_is_reconstructed_when_only_s_and_d_are_reported():
    with_c = dict(PT["lja"])
    without_c = {k: v for k, v in with_c.items() if k != "busco_c"}
    a, _ = score_assembly(with_c, mode="full", expected_haploid=G)
    b, _ = score_assembly(without_c, mode="full", expected_haploid=G)
    assert abs(a - b) < 1000


def test_contributions_sum_to_the_score():
    for mode in DELIVERABLE_MODES:
        total, c = score_assembly(PT["lja"], mode=mode, taxon="fungal",
                                  expected_haploid=G)
        assert abs(total - sum(v for k, v in c.items() if k != "_meta")) < 1e-6


# ── size: reported in full mode, never ranked ────────────────────────────────

def test_size_is_reported_but_not_scored_in_full_mode():
    _, c = score_assembly(PT["lja"], mode="full", taxon="fungal", expected_haploid=G)
    assert c["size_penalty"] == 0
    assert c["_meta"]["size_ranks"] is False
    assert c["_meta"]["size_band_high"] is not None
    _, c = score_assembly(PT["lja"], mode="primary", taxon="fungal", expected_haploid=G)
    assert c["size_penalty"] < 0 and c["_meta"]["size_ranks"] is True


def test_the_artifact_is_not_caught_by_the_size_band_alone():
    """287 Mb sits inside a 2x band on 127 Mb, which is why size cannot rank."""
    deviation, low, high = size_expectation(
        PT["mbg"]["total_len"], G, mode="full")
    assert deviation == 0.0 and low <= PT["mbg"]["total_len"] <= high


def test_size_band_still_penalises_an_assembly_below_the_haploid_size():
    deviation, low, high = size_expectation(60_000_000, G, mode="full")
    assert deviation > 0 and low == G and high > low


def test_contig_density_is_ploidy_neutral():
    _, cl = score_assembly(PT["lja"], mode="full", taxon="fungal", expected_haploid=G)
    _, ci = score_assembly(PT["ipa"], mode="full", taxon="fungal", expected_haploid=G)
    assert cl["contig_count_penalty"] == 0 and ci["contig_count_penalty"] == 0
    assert round(cl["_meta"]["contigs_per_mb"], 1) == 3.8
    assert round(ci["_meta"]["contigs_per_mb"], 1) == 2.5


def test_density_penalty_falls_back_to_absolute_count_without_a_genome_size():
    _, c = score_assembly(PT["lja"], mode="full", taxon="fungal", expected_haploid=0)
    rate = resolve_weights("full", "fungal")["w_contig_density"]
    assert abs(c["contig_density_penalty"] + PT["lja"]["contigs"] * rate) < 1e-6


# ── the telomere reward is bounded by biology ────────────────────────────────

def test_cap_scales_with_ploidy_and_genome_size():
    assert telomere_reward_cap(G, mode="primary") == 254.0
    assert telomere_reward_cap(G, mode="full") == 508.0
    assert telomere_reward_cap(3_000_000_000, mode="primary") == 6000.0


def test_no_cap_without_a_declared_genome_size():
    assert telomere_reward_cap(0) is None
    _, c = score_assembly(PT["mbg"], mode="primary", taxon="fungal", expected_haploid=0)
    assert c["_meta"]["telomere_capped"] is False


def test_the_impossible_telomere_count_is_capped():
    """1,668 telomere-bearing contigs cannot exist in an 18-chromosome genome:
    a contig has two ends, so at most 36 are possible (72 as a dikaryon)."""
    for mode in DELIVERABLE_MODES:
        _, c = score_assembly(PT["mbg"], mode=mode, taxon="fungal", expected_haploid=G)
        assert c["_meta"]["telomere_capped"] is True
        cap = c["_meta"]["telomere_cap"]
        assert c["single_telomere"] == cap * resolve_weights(mode, "fungal")["w_single_tel"]


def test_real_candidates_are_never_capped():
    """The bound is a safety net; it must not touch a plausible assembly."""
    for name, m in PT.items():
        if name == "mbg":
            continue
        for mode in DELIVERABLE_MODES:
            _, c = score_assembly(m, mode=mode, taxon="fungal", expected_haploid=G)
            assert c["_meta"]["telomere_capped"] is False, (name, mode)


# ── what may be removed ──────────────────────────────────────────────────────

def test_primary_permits_every_collapsing_step():
    p = AssemblyPolicy(mode="primary")
    assert p.purge_dups_enabled and p.collapse_redundancy_enabled


def test_full_declines_every_collapsing_step():
    p = AssemblyPolicy(mode="full")
    assert not p.purge_dups_enabled and not p.collapse_redundancy_enabled


def test_contamination_and_chimera_detection_stay_on_in_every_mode():
    for mode in VALID_MODES:
        p = AssemblyPolicy(mode=mode)
        assert p.contamination_screen_enabled and p.chimera_resolution_enabled


def test_no_purge_dups_can_only_disable_never_enable():
    assert not AssemblyPolicy(mode="primary", user_no_purge_dups=True).purge_dups_enabled
    assert not AssemblyPolicy(mode="full", user_no_purge_dups=True).purge_dups_enabled
    assert not AssemblyPolicy(mode="full", user_no_purge_dups=False).purge_dups_enabled
    assert AssemblyPolicy(mode="primary", user_no_purge_dups=False).purge_dups_enabled


def test_no_purge_dups_propagates_into_both_sub_policies():
    for sub in AssemblyPolicy(mode=MODE_BOTH, user_no_purge_dups=True).sub_policies():
        assert not sub.purge_dups_enabled


def test_tool_checks_ask_every_deliverable_not_the_aggregate():
    """check_requirements() must not ask a 'both' policy a removal question.

    Regression: it did, and a default-mode run aborted before step 0 because the
    guard correctly refused to guess. purge_dups is needed when ANY deliverable
    wants it.
    """
    p = AssemblyPolicy(mode=MODE_BOTH)
    assert any(q.purge_dups_enabled for q in p.sub_policies())
    off = AssemblyPolicy(mode=MODE_BOTH, user_no_purge_dups=True)
    assert not any(q.purge_dups_enabled for q in off.sub_policies())
    src = open(os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), "taco", "pipeline.py")).read()
    assert "self.policy.purge_dups_enabled" not in src, (
        "pipeline.py asks the aggregate policy a per-deliverable question")


# ── the combined root report ─────────────────────────────────────────────────

def _fake_both_run(tmp):
    """A finished both-mode run: step-10 table, two deliverables, stale root."""
    import csv as _csv
    os.makedirs(os.path.join(tmp, "assemblies"), exist_ok=True)
    os.makedirs(os.path.join(tmp, "final_results"), exist_ok=True)
    cols = ["ipa", "peregrine", "lja"]
    rows = [("BUSCO C (%)", [96.3, 95.6, 96.7]),
            ("# contigs", [329, 99, 972]),
            ("Total length", [133691058, 147196254, 253991215])]
    with open(os.path.join(tmp, "assemblies", "assembly_info.csv"), "w",
              newline="") as f:
        w = _csv.writer(f)
        w.writerow(["Metric"] + cols)
        for m, v in rows:
            w.writerow([m] + v)
    # a stale root result from an EARLIER run
    with open(os.path.join(tmp, "final_results", "final_result.csv"), "w") as f:
        f.write("Metric,merged\nBUSCO C (%),0.0\n")
    delivered = []
    for mode, pick in (("primary", "peregrine"), ("full", "lja")):
        ws = os.path.join(tmp, f"mode_{mode}")
        os.makedirs(os.path.join(ws, "final_results"), exist_ok=True)
        with open(os.path.join(ws, "final_results", "final_result.csv"), "w",
                  newline="") as f:
            w = _csv.writer(f)
            w.writerow(["Metric", "merged"])
            w.writerow(["BUSCO C (%)", 96.5])
            w.writerow(["# contigs", 104 if mode == "primary" else 940])
            w.writerow(["Selected assembler", pick])
        delivered.append((mode, ws))
    return delivered


def _runner_for_report():
    from taco.pipeline import PipelineRunner

    class R(PipelineRunner):
        def __init__(self):
            self.assembly_mode = MODE_BOTH
            self.taxon = "fungal"
            self.no_purge_dups = False
            self.policy = AssemblyPolicy(mode=MODE_BOTH, taxon="fungal")
            self.run_id = "TEST"
            self.logged = []
            self.benchmark = False
            self.steps = []
            self.logs_dir = "logs"
            self.assembly_only = False
        def log(self, m): self.logged.append(m)
        def log_info(self, m): self.logged.append(m)
        def log_warn(self, m): self.logged.append(m)
        def log_error(self, m): self.logged.append(m)
        # The step loop records timings and restores resume inputs; neither is
        # under test here.
        def append_step_benchmark(self, *a, **k): pass
        def write_benchmark_summary(self): pass
        def restore_resume_inputs_for_step(self, step): pass
        def warn_missing_step_inputs(self, step): pass
    return R()


def test_root_report_carries_every_assembler_and_both_deliverables():
    """One table: the assemblers from step 10 plus a column per deliverable.

    Regression: under 'both' every report was written inside mode_*/ and the
    root final_result.csv was left holding a PREVIOUS run's numbers, which read
    as current for as long as the run took.
    """
    import tempfile
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as tmp:
        delivered = _fake_both_run(tmp)
        try:
            os.chdir(tmp)
            _runner_for_report()._write_both_modes_report(delivered)
            with open(os.path.join("final_results", "final_result.csv"),
                      newline="") as f:
                import csv as _csv
                rows = list(_csv.reader(f))
        finally:
            os.chdir(cwd)
    hdr = rows[0]
    assert hdr == ["Metric", "ipa", "peregrine", "lja",
                   "merged_primary", "merged_full"], hdr
    body = {r[0]: r for r in rows[1:]}
    # the stale value (0.0) must be gone, replaced by both deliverables
    assert body["BUSCO C (%)"][-2:] == ["96.5", "96.5"]
    assert body["BUSCO C (%)"][1] == "96.3"          # assembler column preserved
    # a metric only a deliverable reports still appears
    assert body["Selected assembler"][-2:] == ["peregrine", "lja"]


def test_root_report_removes_the_in_flight_marker():
    import tempfile
    from taco.pipeline import PipelineRunner
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as tmp:
        delivered = _fake_both_run(tmp)
        try:
            os.chdir(tmp)
            r = _runner_for_report()
            r._mark_root_results_in_flight(list(DELIVERABLE_MODES))
            marker = os.path.join("final_results",
                                  PipelineRunner.IN_FLIGHT_MARKER)
            assert os.path.isfile(marker)
            text = open(marker).read()
            # it must name the stale file and where this run is reporting
            assert "final_result.csv" in text
            assert "mode_primary/final_results/" in text
            assert "STILL IN PROGRESS" in text
            r._write_both_modes_report(delivered)
            assert not os.path.exists(marker)
        finally:
            os.chdir(cwd)


def test_deliverable_products_are_never_shared_between_modes():
    """Regression: a real run died on this, and would have corrupted results.

    The mode workspace links step 0-10 products. It also linked busco/final,
    merqury/final and assemblies/final.telo.fasta, which belong to whichever
    representation PRODUCED them. Two consequences: step 13 tried to refresh
    busco/final and hit `Cannot call rmtree on a symbolic link`, and had it not
    crashed, both deliverables would have written their final QC through the same
    links -- silently reporting one representation's metrics for the other.
    """
    from taco.pipeline import PipelineRunner as P
    for name in ("final", "final.telo.fasta", "final.busco.stdout.log",
                 "final_result.csv", "busco/final", "merqury/final"):
        assert P._is_deliverable_product(name), name
    for name in ("canu", "ipa.result.fasta", "compare", "peregrine.telo.fasta",
                 "assembly_info.csv", "finalize.txt", "semifinal"):
        assert not P._is_deliverable_product(name), name


def test_mode_workspace_skips_deliverable_products():
    import tempfile
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as tmp:
        os.makedirs(os.path.join(tmp, "assemblies"), exist_ok=True)
        os.makedirs(os.path.join(tmp, "busco", "final"), exist_ok=True)
        os.makedirs(os.path.join(tmp, "busco", "ipa"), exist_ok=True)
        for n in ("ipa.result.fasta", "ipa.telo.fasta", "final.telo.fasta"):
            open(os.path.join(tmp, "assemblies", n), "w").write(">x\nA\n")
        try:
            os.chdir(tmp)
            r = _runner_for_report()
            ws = r._prepare_mode_workspace("primary")
            assert os.path.lexists(os.path.join(ws, "assemblies", "ipa.result.fasta"))
            assert os.path.lexists(os.path.join(ws, "busco", "ipa"))
            # the deliverable products must be absent, so this run builds its own
            assert not os.path.lexists(os.path.join(ws, "assemblies", "final.telo.fasta"))
            assert not os.path.lexists(os.path.join(ws, "busco", "final"))
        finally:
            os.chdir(cwd)


def test_stale_busco_dir_can_be_a_symlink():
    """The wipe path must not raise when the stale result is a link."""
    src = open(os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), "taco", "steps.py")).read()
    assert "Cannot call rmtree" not in src
    assert src.count("if os.path.islink(") >= 2, (
        "rmtree sites are not guarded against symlinks")


# ── external commands must be time-boundable ─────────────────────────────────

def test_run_cmd_kills_the_whole_process_group_on_timeout():
    """Regression: an optional add-on hung a run and orphaned CPU-burning children.

    dnadiff's delta-filter stage ran 19 h on a 959 MB delta with a zero-byte
    output, blocking step 14 behind an OPTIONAL table; two orphans from earlier
    runs had been going 8 and 11 days. Killing only the immediate child is not
    enough -- a shell pipeline's real worker is a grandchild -- so the whole
    process group must go.
    """
    import subprocess as sp
    import tempfile
    from taco.pipeline import PipelineRunner

    class R(PipelineRunner):
        def __init__(self):
            self.logged = []
        def log(self, m): self.logged.append(m)
        def log_info(self, m): self.logged.append(m)
        def log_warn(self, m): self.logged.append(m)
        def log_error(self, m): self.logged.append(m)

    r = R()
    with tempfile.TemporaryDirectory() as tmp:
        marker = os.path.join(tmp, "still_running")
        # A shell that backgrounds a grandchild which would outlive a naive kill.
        cmd = (f"sh -c 'while :; do echo x >> {marker}; sleep 0.2; done' & "
               f"wait")
        res = r.run_cmd(cmd, check=False, timeout=2)
        assert res.returncode == 124, res.returncode
        assert any("time limit" in m for m in r.logged), r.logged
        size_after_kill = os.path.getsize(marker) if os.path.exists(marker) else 0
        # give any survivor a chance to keep writing
        sp.run(["sh", "-c", "sleep 2"], check=False)
        grew = ((os.path.getsize(marker) if os.path.exists(marker) else 0)
                > size_after_kill)
        assert not grew, "a grandchild survived the timeout kill"


def test_run_cmd_without_a_timeout_is_unchanged():
    from taco.pipeline import PipelineRunner

    class R(PipelineRunner):
        def __init__(self): pass
        def log(self, m): pass
        def log_info(self, m): pass
        def log_warn(self, m): pass
        def log_error(self, m): pass

    assert R().run_cmd("true", check=False).returncode == 0
    assert R().run_cmd("false", check=False).returncode != 0


def test_dnadiff_is_time_bounded_and_optional():
    """It is a convenience table and must never hold up a run."""
    from taco.steps import DNADIFF_TIMEOUT_SEC
    assert 0 < DNADIFF_TIMEOUT_SEC <= 86400
    src = open(os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), "taco", "steps.py")).read()
    call = src[src.index("Optional dnadiff"):]
    call = call[:call.index("def ", 1)] if "def " in call[1:] else call[:4000]
    assert "timeout=DNADIFF_TIMEOUT_SEC" in call
    assert "check=False" in call, "a timed-out optional add-on must not fail the step"


def test_one_failing_deliverable_does_not_cancel_the_other():
    """Regression: it did, and it cost a day of compute.

    A step-13 crash in the primary deliverable aborted the whole run, so
    mode_full never started. The representations are independent products and a
    run of this length must salvage the one that works. _run_step_list exits
    rather than raising on a failed step, and SystemExit does not derive from
    Exception, so it has to be caught explicitly.
    """
    import tempfile
    cwd = os.getcwd()
    attempted = []

    def primary_explodes(runner):
        attempted.append(runner.assembly_mode)
        if runner.assembly_mode == "primary":
            raise SystemExit(1)            # exactly how a failed step reports
        os.makedirs("final_results", exist_ok=True)
        import csv as _csv
        with open(os.path.join("final_results", "final_result.csv"), "w",
                  newline="") as f:
            w = _csv.writer(f)
            w.writerow(["Metric", "merged"])
            w.writerow(["Selected assembler", "lja"])

    with tempfile.TemporaryDirectory() as tmp:
        os.makedirs(os.path.join(tmp, "assemblies"))
        open(os.path.join(tmp, "assemblies", "ipa.result.fasta"), "w").write(">x\nA\n")
        r = _runner_for_report()
        r.steps = [11]
        r.logs_dir = os.path.join(tmp, "logs")
        try:
            os.chdir(tmp)
            try:
                r._run_both_modes({11: primary_explodes})
            except SystemExit as e:
                # a failed deliverable must still make the RUN fail...
                assert e.code == 1, e.code
            else:
                raise AssertionError("a failed deliverable must exit non-zero")
            # ...but only after the other one was attempted and reported
            assert attempted == ["primary", "full"], attempted
            out = os.path.join("final_results", "final_result.csv")
            assert os.path.isfile(out), "the surviving deliverable was not reported"
            text = open(out).read()
            assert "merged_full" in text, text
            assert "merged_primary" not in text, text
        finally:
            os.chdir(cwd)


def test_timeout_kill_reads_the_process_group_once():
    """Regression: re-reading the pgid after signalling raced with process exit,
    turning a handled timeout into an uncaught ProcessLookupError."""
    src = open(os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), "taco", "pipeline.py")).read()
    body = src[src.index("def _kill_process_group"):]
    body = body[:body.index("\n    def ", 1)]
    assert body.count("os.getpgid(") == 1, (
        "the process group id must be read once, before signalling")
    assert "os.killpg(pgid" in body


def test_run_cmd_accepts_a_list_command_with_a_timeout():
    from taco.pipeline import PipelineRunner

    class R(PipelineRunner):
        def __init__(self): pass
        def log(self, m): pass
        def log_info(self, m): pass
        def log_warn(self, m): pass
        def log_error(self, m): pass

    assert R().run_cmd(["true"], check=False).returncode == 0
    assert R().run_cmd(["sleep", "10"], check=False, timeout=1).returncode == 124


def test_combined_report_survives_a_missing_deliverable_csv():
    """A partially completed run must still produce a readable root table."""
    import tempfile
    import csv as _csv
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as tmp:
        os.makedirs(os.path.join(tmp, "assemblies"))
        os.makedirs(os.path.join(tmp, "mode_primary", "final_results"))
        with open(os.path.join(tmp, "assemblies", "assembly_info.csv"), "w",
                  newline="") as f:
            w = _csv.writer(f)
            w.writerow(["Metric", "ipa"])
            w.writerow(["BUSCO C (%)", 96.3])
        with open(os.path.join(tmp, "mode_primary", "final_results",
                               "final_result.csv"), "w", newline="") as f:
            w = _csv.writer(f)
            w.writerow(["Metric", "merged"])
            w.writerow(["BUSCO C (%)", 92.5])
        try:
            os.chdir(tmp)
            _runner_for_report()._write_both_modes_report(
                [("primary", os.path.join(tmp, "mode_primary")),
                 ("full", os.path.join(tmp, "mode_full"))])   # full never ran
            rows = list(_csv.reader(open(os.path.join("final_results",
                                                      "final_result.csv"))))
        finally:
            os.chdir(cwd)
    assert rows[0][-2:] == ["merged_primary", "merged_full"], rows[0]
    assert rows[1][-2:] == ["92.5", ""], rows[1]


def test_mode_workspace_preparation_is_idempotent():
    """Resuming a deliverable re-runs it; it must not fail on existing links."""
    import tempfile
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as tmp:
        os.makedirs(os.path.join(tmp, "assemblies"))
        open(os.path.join(tmp, "assemblies", "ipa.result.fasta"), "w").write(">x\nA\n")
        try:
            os.chdir(tmp)
            r = _runner_for_report()
            a = r._prepare_mode_workspace("full")
            b = r._prepare_mode_workspace("full")
            assert a == b
            assert os.path.lexists(os.path.join(b, "assemblies", "ipa.result.fasta"))
        finally:
            os.chdir(cwd)


# ── purge_dups must not remove complete genes unnoticed ──────────────────────
#
# Real Puccinia triticina data: purging a peregrine backbone removed 39
# haplotigs and 25.7 Mb, took BUSCO complete from 95.6% to 92.5%, and passed
# every existing check. The 12L do-no-harm test compares assembly size and
# telomere counts -- and its size rule explicitly APPROVES a shrink that moves
# toward the declared genome size, which is exactly what a purge does.

def test_busco_thresholds_are_shared_by_every_harm_check():
    from taco.steps import _busco_delta_thresholds
    assert _busco_delta_thresholds("fungal")[0] == "2.0"
    assert _busco_delta_thresholds("plant")[0] == "4.0"
    assert _busco_delta_thresholds("other") == _busco_delta_thresholds("insect")
    src = open(os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), "taco", "steps.py")).read()
    # The table must exist in exactly one place.
    assert src.count('return "4.0", "1.0", "6.0"') == 1
    assert src.count("_busco_delta_thresholds(") >= 3, "gate/trial/report share it"


def test_busco_complete_is_read_from_the_summary_busco_already_wrote():
    import tempfile
    from taco.steps import _busco_complete_for_label
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as tmp:
        d = os.path.join(tmp, "busco", "final")
        os.makedirs(d)
        with open(os.path.join(d, "short_summary.specific.x.final.txt"), "w") as f:
            f.write("\tC:92.5%[S:87.4%,D:5.1%],F:0.6%,M:6.9%,n:1752\n")
        try:
            os.chdir(tmp)
            assert _busco_complete_for_label("final") == 92.5
            assert _busco_complete_for_label("absent") is None
        finally:
            os.chdir(cwd)


def _gate_runner(taxon="fungal"):
    class R:
        def __init__(self):
            self.taxon = taxon
            self.busco_lineage = "basidiomycota_odb10"
            self.threads = 4
            self.logged = []
        def log(self, m): self.logged.append(m)
        def log_info(self, m): self.logged.append(m)
        def log_warn(self, m): self.logged.append(m)
    return R()


def _run_gate(before_c, after_c, taxon="fungal", before_d=20.0, after_d=5.0):
    """Drive the gate with mocked BUSCO trials."""
    from taco import steps as st
    import tempfile
    real_trial, real_which = st._run_busco_trial, st.shutil.which
    seq = {"purge_before": {"C": before_c, "D": before_d},
           "purge_after": {"C": after_c, "D": after_d}}
    st._run_busco_trial = lambda fa, lin, thr, label, out, runner=None: seq[label]
    st.shutil.which = lambda n: "/usr/bin/busco" if n == "busco" else real_which(n)
    cwd = os.getcwd()
    try:
        with tempfile.TemporaryDirectory() as tmp:
            os.chdir(tmp)
            r = _gate_runner(taxon)
            return st._purge_preserves_gene_content(r, "before.fa", "after.fa"), r
    finally:
        os.chdir(cwd)
        st._run_busco_trial, st.shutil.which = real_trial, real_which


def test_gate_rejects_the_observed_gene_loss():
    """95.6 -> 92.5 is a 3.1-point drop against a 2.0-point fungal tolerance."""
    ok, r = _run_gate(95.6, 92.5)
    assert ok is False
    assert any("92.5" in m for m in r.logged), r.logged


def test_gate_accepts_a_purge_that_only_removes_redundancy():
    ok, r = _run_gate(96.3, 96.2)
    assert ok is True
    assert any("accepted" in m for m in r.logged), r.logged


def test_gate_tolerance_follows_the_taxon():
    """The same 3.1-point drop is inside a plant tolerance of 4.0."""
    assert _run_gate(95.6, 92.5, taxon="plant")[0] is True
    assert _run_gate(95.6, 92.5, taxon="fungal")[0] is False


def test_gate_can_be_overridden_but_says_so():
    os.environ["STEP12_SKIP_PURGE_BUSCO_GATE"] = "1"
    try:
        ok, r = _run_gate(95.6, 50.0)
        assert ok is True
        assert any("gate skipped" in m for m in r.logged), r.logged
    finally:
        del os.environ["STEP12_SKIP_PURGE_BUSCO_GATE"]


def test_gate_is_permissive_when_busco_cannot_measure():
    """Missing BUSCO must not silently block a purge, but must warn."""
    from taco import steps as st
    real = st._run_busco_trial
    st._run_busco_trial = lambda *a, **k: None
    try:
        r = _gate_runner()
        assert st._purge_preserves_gene_content(r, "b.fa", "a.fa") is True
        assert any("inconclusive" in m or "unavailable" in m for m in r.logged)
    finally:
        st._run_busco_trial = real


def test_a_rejected_purge_keeps_the_unpurged_assembly():
    src = open(os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), "taco", "steps.py")).read()
    blk = src[src.index("if _pol.purge_dups_enabled:"):]
    blk = blk[:blk.index("elif _pol.user_no_purge_dups")]
    assert "_purge_preserves_gene_content" in blk
    # the copy that overwrites the assembly must be conditional on the gate
    assert blk.index("_purge_preserves_gene_content") < blk.index("shutil.copy")
    assert "REJECTED" in blk


def test_do_no_harm_checks_gene_content():
    src = open(os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), "taco", "steps.py")).read()
    blk = src[src.index('# ---- 12L. "Do no harm" safety comparison ----'):]
    blk = blk[:blk.index("\ndef ", 1)]
    assert "_busco_complete_for_label" in blk, "12L still ignores gene content"
    assert "BUSCO complete" in blk


# ── full mode needs candidates that actually retain alternate sequence ───────
#
# Several assemblers emit more than one contig set and TACO used only the
# primary, so `--assembly-mode full` chose among assemblies that had all been
# reduced to one representation before it saw them. On real P. triticina data
# hifiasm and IPA reported 5.2% and 6.9% BUSCO duplication while canu and LJA
# reported 91.2% and 95.9% -- the full-representation winner was decided by
# which assemblers collapse by default, not by which recovers alternate
# sequence best. hifiasm had already written bp.hap1.p_ctg and bp.hap2.p_ctg,
# and IPA a_ctg.fasta; both were discarded.

def test_derived_candidates_are_registered_with_their_parent():
    from taco.utils import ALL_ASSEMBLERS, DERIVED_ASSEMBLERS
    assert DERIVED_ASSEMBLERS == {"hifiasm_hap": "hifiasm", "ipa_alt": "ipa"}
    for child, parent in DERIVED_ASSEMBLERS.items():
        assert child in ALL_ASSEMBLERS, child
        assert parent in ALL_ASSEMBLERS, parent


def test_a_derived_candidate_does_not_vote_beside_its_parent():
    """One vote per INDEPENDENT assembly: a derived set shares its parent's
    graph, reads and error profile, so two ballots would inflate agreement."""
    from taco.utils import concordance_voters
    got = concordance_voters(["canu", "hifiasm", "hifiasm_hap", "ipa",
                              "ipa_alt", "lja"])
    assert got == ["canu", "hifiasm", "ipa", "lja"], got
    # parents keep voting even when the child is absent
    assert concordance_voters(["hifiasm", "canu"]) == ["hifiasm", "canu"]


def test_derived_candidates_are_still_selectable():
    """They must compete for the backbone -- only the VOTE is withheld."""
    from taco.utils import DERIVED_ASSEMBLERS, EXCLUDED_FROM_BACKBONE
    for child in DERIVED_ASSEMBLERS:
        assert child not in EXCLUDED_FROM_BACKBONE, child


def test_concordance_applies_the_voter_filter():
    src = open(os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), "taco", "steps.py")).read()
    blk = src[src.index("def _t2t_concordance_check"):]
    blk = blk[:blk.index("\ndef ", 1)]
    assert "asms = concordance_voters(asms)" in blk, (
        "the voter filter is imported but never applied")


def test_normalize_knows_where_the_derived_outputs_live():
    src = open(os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), "taco", "steps.py")).read()
    for path in ("hifiasm/hifiasm.hap.fasta",
                 "ipa/assembly-results/ipa.alt.fasta"):
        assert path in src, path


def test_hifiasm_haplotype_set_is_built_from_both_haplotypes():
    """It must combine hap1 AND hap2; either alone is not the full set."""
    src = open(os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), "taco", "steps.py")).read()
    fn = src[src.index("def _emit_hifiasm_haplotype_set"):]
    fn = fn[:fn.index("\ndef ", 1)]
    assert 'for hap in ("hap1", "hap2")' in fn
    assert "bp.{hap}.p_ctg.gfa" in fn


def test_concat_fasta_never_fuses_records_across_a_missing_newline():
    """A source FASTA without a trailing newline would otherwise glue its last
    sequence line onto the next file's header, producing one fused contig with a
    '>' buried in its sequence. A contig count alone would not reveal it."""
    import tempfile
    from taco import steps as st

    class R:
        def log_info(self, m): pass
        def log_warn(self, m): pass

    cases = [(">h1tg1\nACGTACGT",   ">h2tg1\nTTTTGGGG\n"),   # first lacks one
             (">h1tg1\nACGTACGT\n", ">h2tg1\nTTTTGGGG\n"),   # both fine
             (">h1tg1\nACGTACGT",   ">h2tg1\nTTTTGGGG")]      # neither
    for a_txt, b_txt in cases:
        with tempfile.TemporaryDirectory() as tmp:
            a = os.path.join(tmp, "a.fa")
            b = os.path.join(tmp, "b.fa")
            out = os.path.join(tmp, "out.fa")
            open(a, "w").write(a_txt)
            open(b, "w").write(b_txt)
            assert st._concat_fasta(R(), [a, b], out, "t") is True
            text = open(out).read()
            lines = text.split("\n")
            assert len([l for l in lines if l.startswith(">")]) == 2, text
            fused = [l for l in lines if ">" in l and not l.startswith(">")]
            assert not fused, f"header fused into sequence: {fused}"


def test_concat_fasta_combines_only_what_exists():
    import tempfile
    from taco import steps as st

    class R:
        def log_info(self, m): pass
        def log_warn(self, m): pass

    with tempfile.TemporaryDirectory() as tmp:
        a = os.path.join(tmp, "a.fa")
        b = os.path.join(tmp, "b.fa")
        out = os.path.join(tmp, "out.fa")
        open(a, "w").write(">h1tg1\nACGT\n")
        # b deliberately absent
        assert st._concat_fasta(R(), [a, b], out, "t") is True
        assert open(out).read() == ">h1tg1\nACGT\n"
        # nothing present at all -> no file claimed
        assert st._concat_fasta(R(), [b], os.path.join(tmp, "x.fa"), "t") is False


# ── the merge is where gene content actually goes ────────────────────────────
#
# Root cause of a 6.3-point BUSCO drop on real Puccinia triticina data. Step 12
# rebuilds the assembly from the telomere pool and drops backbone sequence the
# pool appears to cover. Every gate sat on a different operation: the 12F trial
# gates individual pool REPLACEMENTS (one upgrade, ten additions on that run)
# and the 12H gate covers purge_dups (which full mode skips entirely), while the
# merge turned 972 backbone contigs into 468. Of the 112 complete BUSCOs lost,
# 104 were MISSING and only 8 Fragmented -- sequence gone, not junction-broken.

def test_busco_full_table_parser_handles_duplicated_copies():
    """full_table lists one line per COPY, so a Duplicated gene has >1 contig."""
    import tempfile
    from taco.steps import _busco_complete_by_id
    with tempfile.TemporaryDirectory() as tmp:
        run = os.path.join(tmp, "run_x")
        os.makedirs(run)
        with open(os.path.join(run, "full_table.tsv"), "w") as f:
            f.write("# comment line\n")
            f.write("g1\tComplete\tctg_1\t1\t100\n")
            f.write("g2\tDuplicated\tctg_2\t1\t100\n")
            f.write("g2\tDuplicated\tctg_3\t1\t100\n")
            f.write("g3\tFragmented\tctg_4\t1\t100\n")
            f.write("g4\tMissing\t\t\t\n")
        got = _busco_complete_by_id(tmp)
    assert got == {"g1": {"ctg_1"}, "g2": {"ctg_2", "ctg_3"}}, got
    assert "g3" not in got and "g4" not in got


def _recovery_runner(taxon="fungal"):
    class R:
        def __init__(self):
            self.taxon = taxon
            self.busco_lineage = "basidiomycota_odb10"
            self.threads = 2
            self.logged = []
        def log(self, m): self.logged.append(m)
        def log_info(self, m): self.logged.append(m)
        def log_warn(self, m): self.logged.append(m)
    return R()


def _drive_recovery(c_backbone, c_merged, backbone_extra_contig=True):
    """Run the recovery pass with mocked BUSCO, return (runner, merged text)."""
    import tempfile
    from taco import steps as st
    real_trial, real_ids, real_which = (
        st._run_busco_trial, st._busco_complete_by_id, st.shutil.which)
    st.shutil.which = lambda n: "/usr/bin/busco" if n == "busco" else real_which(n)
    st._run_busco_trial = lambda fa, lin, thr, label, out, runner=None: (
        {"C": c_backbone} if label == "merge_backbone" else {"C": c_merged})
    st._busco_complete_by_id = lambda d: (
        {"g1": {"bb_1"}, "g2": {"bb_2"}} if d.endswith("merge_backbone")
        else {"g1": {"mg_1"}})
    cwd = os.getcwd()
    try:
        tmp = tempfile.mkdtemp()
        os.chdir(tmp)
        os.makedirs("assemblies", exist_ok=True)
        bb, mg = "assemblies/bb.fa", "assemblies/merged.fa"
        with open(bb, "w") as f:
            f.write(">bb_1\nACGT\n")
            if backbone_extra_contig:
                f.write(">bb_2\nTTTTGGGG\n")
        open(mg, "w").write(">mg_1\nACGT\n")
        r = _recovery_runner()
        st._recover_genes_lost_in_merge(r, bb, mg)
        return r, open(mg).read()
    finally:
        os.chdir(cwd)
        st._run_busco_trial, st._busco_complete_by_id, st.shutil.which = (
            real_trial, real_ids, real_which)


def test_merge_that_loses_genes_gets_them_restored():
    """96.7 -> 90.4 is 6.3 points, past the 2.0-point fungal tolerance."""
    r, merged = _drive_recovery(96.7, 90.4)
    assert "bb_2_gene_rescue" in merged, merged
    assert "TTTTGGGG" in merged
    assert any("Restored" in m for m in r.logged), r.logged


def test_a_merge_that_preserves_genes_restores_nothing():
    r, merged = _drive_recovery(96.7, 96.6)      # 0.1 drop, inside tolerance
    assert "gene_rescue" not in merged, merged
    assert any("nothing to recover" in m for m in r.logged), r.logged


def test_recovery_is_reported_even_when_no_contig_can_be_identified():
    """A loss with no dropped carrier must warn rather than fail silently."""
    r, merged = _drive_recovery(96.7, 90.4, backbone_extra_contig=False)
    assert "gene_rescue" not in merged
    assert any("not explained by dropped contigs" in m for m in r.logged), r.logged


def test_recovery_can_be_disabled_but_says_so():
    os.environ["STEP12_SKIP_MERGE_GENE_RECOVERY"] = "1"
    try:
        r, merged = _drive_recovery(96.7, 50.0)
        assert "gene_rescue" not in merged
        assert any("skipped" in m for m in r.logged), r.logged
    finally:
        del os.environ["STEP12_SKIP_MERGE_GENE_RECOVERY"]


def test_the_recovery_pass_runs_before_purge_dups():
    """It must see the merge output, not a purged version of it."""
    src = open(os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), "taco", "steps.py")).read()
    assert src.index("_recover_genes_lost_in_merge(") < src.index(
        "if _pol.purge_dups_enabled:")


def test_every_content_removing_stage_now_has_a_gene_gate():
    """Merge, purge_dups, and the final comparison each check gene content."""
    src = open(os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), "taco", "steps.py")).read()
    for guard in ("_recover_genes_lost_in_merge",      # 12G3, the merge
                  "_purge_preserves_gene_content",     # 12H, purge_dups
                  "_busco_complete_for_label"):        # 12L, the report
        assert guard in src, guard


def test_policy_rejects_an_unknown_mode():
    try:
        AssemblyPolicy(mode="phased")
    except ValueError:
        pass
    else:
        raise AssertionError("AssemblyPolicy accepted an unknown mode")


def test_full_mode_tightens_the_shallow_depth_threshold():
    """A retained haplotype sits near 0.5x core depth; 0.34 leaves only 1.5x."""
    from taco.purify import DEPTH_RATIO_MAX, DEPTH_RATIO_MAX_FULL
    assert DEPTH_RATIO_MAX_FULL < DEPTH_RATIO_MAX
    assert DEPTH_RATIO_MAX_FULL < 0.5 / 1.9


def test_no_mode_is_ever_described_as_phased():
    """The feature retains alternate sequence; it does not assign haplotypes."""
    banned = ("hap1", "hap2", "haplotype-resolved", "haplotype resolved")
    surfaces = [AssemblyPolicy(mode=m).describe() for m in VALID_MODES]
    surfaces += [SCORING_PROFILES[m]["description"] for m in DELIVERABLE_MODES]
    for text in surfaces:
        low = text.lower()
        for term in banned:
            assert term not in low, f"{term!r} used to describe a mode: {text!r}"
        assert "phased" not in low, text


# ── advisory ─────────────────────────────────────────────────────────────────

def test_advisory_fires_on_the_three_condition_pattern():
    msg = haplotype_retention_advisory(PT["lja"], G, mode="primary")
    assert msg is not None
    assert "do not establish the ploidy" in msg


def test_advisory_makes_no_ploidy_claim():
    msg = haplotype_retention_advisory(PT["lja"], G, mode="primary").lower()
    for term in ("is diploid", "is dikaryotic", "confirms", "proves"):
        assert term not in msg, term
    assert "may indicate" in msg


def test_advisory_is_silent_when_a_full_genome_is_already_being_delivered():
    for mode in ("full", MODE_BOTH):
        assert haplotype_retention_advisory(PT["lja"], G, mode=mode) is None


def test_advisory_needs_all_three_signals():
    assert haplotype_retention_advisory(PT["ipa"], G, mode="primary") is None
    for tweak in (dict(busco_d=1.0), dict(merqury_comp=70.0),
                  dict(total_len=130_000_000)):
        assert haplotype_retention_advisory(
            dict(PT["lja"], **tweak), G, mode="primary") is None


def test_advisory_is_silent_without_a_genome_size():
    assert haplotype_retention_advisory(PT["lja"], 0, mode="primary") is None


# ── backbone.py delegates to the same policy ─────────────────────────────────

BB = {name: {"BUSCO_C": m["busco_c"], "BUSCO_S": m["busco_s"], "BUSCO_D": m["busco_d"],
             "MerquryComp": m["merqury_comp"], "MerquryQV": m["merqury_qv"],
             "contigs": m["contigs"], "N50": m["n50"], "T2T": m["t2t"],
             "single": m["single_tel"], "total_length": m["total_len"]}
      for name, m in PT.items()}


def test_backbone_module_agrees_with_the_policy():
    """backbone.py used to carry its own divergent copy of the formula."""
    for name in BB:
        for mode in DELIVERABLE_MODES:
            expected, _ = score_assembly(PT[name], mode=mode, taxon="fungal",
                                         expected_haploid=G)
            got = _compute_score(BB[name], taxon="fungal", genomesize=G,
                                 assembly_mode=mode)
            assert abs(expected - got) < 1e-6, (name, mode)


def test_backbone_selection_follows_the_mode():
    assert select_backbone(BB, taxon="fungal", genomesize=G,
                           assembly_mode="primary") == "peregrine"
    assert select_backbone(BB, taxon="fungal", genomesize=G,
                           assembly_mode="full") == "lja"


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
