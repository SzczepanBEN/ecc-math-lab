"""Tests for the X_1(13) Tate-normal-form model (experiments/phase3)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments" / "phase3"))

from tate13 import F, PTS, pari, reduce_mod_F, relation_matrix  # noqa: E402

b, c = pari("b"), pari("c")


def test_small_multiples_match_known_formulas():
    assert PTS[2] == (b, b * c)          # 2P = (b, bc)
    assert PTS[3] == (c, b - c)          # 3P = (c, b - c)


def test_order_13_locus():
    fa = pari.factor(F)
    assert len(fa[0]) == 1 and fa[1][0] == 1
    assert pari.pollead(F, "b") == 1
    assert (int(pari.poldegree(F, "b")), int(pari.poldegree(F, "c"))) == (7, 10)


def test_T13_core_identity():
    assert reduce_mod_F(pari.numerator(pari.matdet(relation_matrix(4, 2)))) == 0
    assert reduce_mod_F(pari.numerator(pari.matdet(relation_matrix(2, 2)))) != 0
