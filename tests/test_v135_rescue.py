"""Regression tests for the v1.3.5 rescue + purge_dups safety changes.

Run from the repo root:  python tests/test_v135_rescue.py
These use only the standard library (no assemblers/BUSCO/minimap2 required).
"""
import os
import sys
import random
import tempfile

# Make `taco` importable when run directly from a checkout.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from taco import steps  # noqa: E402


def test_overhang_telomere_detector():
    telo = "TTAGGG" * 30
    rnd = "".join(random.Random(1).choice("ACGT") for _ in range(180))
    t_telo, s_telo = steps._overhang_is_telomeric(telo, taxon="fungal")
    t_rnd, _ = steps._overhang_is_telomeric(rnd, taxon="fungal")
    assert t_telo is True and s_telo > 0.08
    assert t_rnd is False


def test_short_telomere_cap_is_rescued():
    """A donor that adds only a short (180 bp) telomere cap must survive the
    minimum-extension gate; a non-telomere short overhang must still be dropped."""
    rng = random.Random(7)
    bb = "".join(rng.choice("ACGT") for _ in range(20000))
    donor_seqs = {"d1": bb + "TTAGGG" * 30,
                  "d2": bb + "".join(rng.choice("ACGT") for _ in range(180))}
    backbone_seqs = {"bb1": bb}
    d = tempfile.mkdtemp()
    paf = os.path.join(d, "r.paf")
    with open(paf, "w") as f:
        f.write("d1\t20180\t0\t20000\t+\tbb1\t20000\t0\t20000\t20000\t20000\t60\n")
        f.write("d2\t20180\t0\t20000\t+\tbb1\t20000\t0\t20000\t20000\t20000\t60\n")
    hits = steps._parse_paf_rescue_hits(paf, donor_seqs, backbone_seqs, taxon="fungal")
    rejected, candidates = steps._screen_rescue_candidates(hits)
    cand = {c["donor"] for c in candidates}
    rej = {r["donor"]: r["reject_reason"] for r in rejected}
    assert "d1" in cand, "telomere-adding short extension was wrongly rejected"
    assert "d2" in rej and "ext=" in rej["d2"], "non-telomere short ext should be rejected"


def test_purge_dups_overshoot_guard():
    class R:
        genomesize = "1000"
        taxon = "fungal"

        def log_warn(self, *a, **k):
            pass

    prof = steps._purge_dups_profile("fungal")
    rep = os.path.join(tempfile.mkdtemp(), "safety.tsv")
    # well-sized (1000) purged down to 880 (further below expected) -> reject
    ok_bad, why_bad = steps._purge_dups_safety_check(
        R(), [("c", "A" * 1000)], [("c", "A" * 880)], prof, prof["mode"], rep)
    # inflated (1200) purged to 1000 (toward expected) -> accept
    ok_ok, _ = steps._purge_dups_safety_check(
        R(), [("c", "A" * 1200)], [("c", "A" * 1000)], prof, prof["mode"], rep)
    assert ok_bad is False and why_bad == "overshoot_below_expected_size"
    assert ok_ok is True


if __name__ == "__main__":
    fails = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"PASS  {name}")
            except AssertionError as e:
                fails += 1
                print(f"FAIL  {name}: {e}")
    print("\n%s" % ("ALL TESTS PASSED" if not fails else f"{fails} TEST(S) FAILED"))
    sys.exit(1 if fails else 0)
