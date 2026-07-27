"""Regression tests for the v1.3.6 correctness review fixes.

Run from the repo root:  python tests/test_review_fixes.py
Uses only the standard library — no assemblers, BUSCO, or minimap2 required,
so this suite is safe to run in CI.
"""
import gzip
import os
import random
import sys
import tempfile

# Make `taco` importable when run directly from a checkout.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from taco import steps  # noqa: E402
from taco import telomere_detect as td  # noqa: E402
from taco.utils import revcomp  # noqa: E402
from taco.pipeline import PipelineRunner  # noqa: E402


def _fams(taxon="fungal"):
    keys = td.get_taxon_families(taxon)
    return {k: td.MOTIF_FAMILIES[k] for k in keys if k in td.MOTIF_FAMILIES}


def test_telomere_end_scoring_is_symmetric():
    """A telomere shorter than the score window must score the same at the 5'
    and 3' terminus. Before v1.3.6 both ends were measured from the left edge
    of the window, which systematically under-scored 3' telomeres."""
    fams = _fams("fungal")
    sw = 300
    rng = random.Random(5)
    core = "".join(rng.choice("ACGT") for _ in range(20000))
    for arr in (90, 150, 300, 480):
        telo = "TTAGGG" * (arr // 6)
        left_seq = td.extract_ends(telo + core, sw)[0]
        right_seq = td.extract_ends(core + telo, sw)[1]
        ls = td.score_end(left_seq, sw, [], fams, end="left")["raw_score"]
        rs = td.score_end(right_seq, sw, [], fams, end="right")["raw_score"]
        assert abs(ls - rs) < 1e-9, f"asymmetry at {arr} bp: L={ls} R={rs}"


def test_short_double_telomere_is_strict_t2t():
    """Short telomere arrays at both ends classify as strict_t2t; one end
    only classifies as single_tel_strong."""
    rng = random.Random(11)
    core = "".join(rng.choice("ACGT") for _ in range(40000))
    telo = "TTAGGG" * 15  # 90 bp, shorter than the 300 bp fungal window
    seqs = {
        "c_t2t": telo + core + telo,
        "c_left": telo + core,
        "c_right": core + telo,
        "c_none": core,
    }
    d = tempfile.mkdtemp()
    fa = os.path.join(d, "t.fa")
    td.write_fasta(seqs, fa)
    res = {r["contig"]: r["classification"]
           for r in td.detect_telomeres(fa, mode="hybrid", end_window=5000,
                                        score_window=300, taxon="fungal")}
    assert res["c_t2t"] == "strict_t2t", res
    assert res["c_left"] == "single_tel_strong", res
    assert res["c_right"] == "single_tel_strong", res
    assert res["c_none"] == "none", res


def test_classify_thresholds_are_applied():
    strong = {"raw_score": 0.30}
    weak = {"raw_score": 0.10}
    none = {"raw_score": 0.01}
    assert td.classify_contig(strong, strong) == "strict_t2t"
    assert td.classify_contig(strong, weak) == "single_tel_strong"
    assert td.classify_contig(weak, none) == "telomere_supported"
    assert td.classify_contig(none, none) == "none"
    # explicit thresholds must be honored
    assert td.classify_contig(weak, weak, strong_thr=0.05, weak_thr=0.02) == "strict_t2t"


def test_revcomp_iupac_and_case():
    assert revcomp("ACGT") == "ACGT"
    assert revcomp("AAAA") == "TTTT"
    assert revcomp("RYKM") == "KMRY"          # R<->Y, K<->M
    assert revcomp("acgtn") == "nacgt"        # lowercase reverse complement
    assert revcomp("ACGT-N") == "N-ACGT"      # non-IUPAC characters preserved


def test_collision_free_name_protects_other_contigs():
    """A donor whose name collides with a DIFFERENT backbone contig must never
    overwrite it — that could silently destroy a protected Tier 1 contig."""
    existing = {"contig_1": "x", "contig_2": "y", "flye_3": "z"}
    # the contig being intentionally replaced keeps its name
    assert steps._collision_free_name("contig_2", existing, exclude="contig_2") == "contig_2"
    # a collision with a different contig is renamed
    got = steps._collision_free_name("contig_1", existing, exclude="flye_3")
    assert got != "contig_1" and got not in existing
    # no collision -> unchanged
    assert steps._collision_free_name("pool_9", existing) == "pool_9"


def test_maybe_gunzip_detects_by_content():
    """A gzip file named .fasta is still decompressed (magic-byte sniff)."""
    d = tempfile.mkdtemp()
    mislabeled = os.path.join(d, "reference_input.fasta")  # gzip bytes, .fasta name
    with gzip.open(mislabeled, "wt") as f:
        f.write(">chr1\nACGTACGTACGT\n")
    runner = PipelineRunner.__new__(PipelineRunner)  # no __init__ needed
    runner.log_info = lambda *a, **k: None
    out = runner._maybe_gunzip(mislabeled)
    assert out != mislabeled
    with open(out) as f:
        assert f.read().startswith(">chr1")
    # a plain-text FASTA is returned unchanged
    plain = os.path.join(d, "plain.fasta")
    with open(plain, "w") as f:
        f.write(">c\nACGT\n")
    assert runner._maybe_gunzip(plain) == plain


def test_busco_full_table_counts_genes_not_copies():
    """A Duplicated BUSCO spans several copy lines in full_table.tsv; counts
    must be per gene so that percentages match BUSCO's own short_summary."""
    lines = ["# Busco id\tStatus\tSequence\tStart\tEnd"]
    for i in range(900):
        lines.append(f"C{i}\tComplete\tctg\t1\t100")
    for i in range(50):  # 50 duplicated genes, 2 copies each = 100 lines
        lines.append(f"D{i}\tDuplicated\tctgA\t1\t100")
        lines.append(f"D{i}\tDuplicated\tctgB\t1\t100")
    for i in range(30):
        lines.append(f"F{i}\tFragmented\tctg\t1\t100")
    for i in range(20):
        lines.append(f"M{i}\tMissing\t\t\t")

    status_idx = 1
    status_by_id = {}
    id_idx = 0 if status_idx != 0 else 1
    for ln in lines:
        if not ln or ln.startswith("#"):
            continue
        parts = ln.split("\t")
        st = parts[status_idx].strip() if len(parts) > status_idx else ""
        bid = parts[id_idx].strip() if len(parts) > id_idx else ""
        if not bid:
            continue
        if st == "Duplicated" or status_by_id.get(bid) == "Duplicated":
            status_by_id[bid] = "Duplicated"
        else:
            status_by_id.setdefault(bid, st)
    n = len(status_by_id)
    D = sum(1 for v in status_by_id.values() if v == "Duplicated")
    S = sum(1 for v in status_by_id.values() if v == "Complete")
    assert n == 1000, f"denominator inflated by copy lines: {n}"
    assert D == 50, f"duplicated counted per copy, not per gene: {D}"
    assert S == 900 and abs(100.0 * (S + D) / n - 95.0) < 1e-9


def test_equal_length_dedup_keeps_one_representative():
    """minimap2 self-alignment emits both A->B and B->A. For an equal-length
    duplicate pair, exactly one member may be dropped — never both."""
    def clean_contained(pairs, pct_cov=30, id_thr=0.90):
        remove = set()
        for q, t, ql, tl, cov, idt in pairs:
            if cov * 100 >= pct_cov and idt >= id_thr and (
                    ql < tl or (ql == tl and q < t)):
                remove.add(q)
        return remove

    def self_dedup(pairs, cov_thr=0.80, id_thr=0.90):
        drop = set()
        for q, t, ql, tl, cov, idt in pairs:
            if cov >= cov_thr and idt >= id_thr:
                if ql < tl or (ql == tl and q < t):
                    drop.add(q)
                else:
                    drop.add(t)
        return drop

    equal = [("A", "B", 1000, 1000, 1.0, 0.99), ("B", "A", 1000, 1000, 1.0, 0.99)]
    assert len(clean_contained(equal)) == 1, "both copies of an equal pair removed"
    assert len(self_dedup(equal)) == 1, "both copies of an equal pair dropped"

    unequal = [("A", "B", 800, 1000, 0.95, 0.98), ("B", "A", 1000, 800, 0.76, 0.98)]
    assert clean_contained(unequal) == {"A"}, "the shorter contig should be removed"

    # a shared repeat at low identity must not trigger removal
    repeat = [("A", "B", 900, 1000, 0.40, 0.75)]
    assert clean_contained(repeat) == set(), "low-identity repeat triggered removal"


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
