"""
Padé Approximant Test
"""
import pytest

from contfrac.pade import PadeApprox


def test_basic():
    "test basic function approximation"
    c = [1,-3/4,39/32]

    p = PadeApprox(1, 1, c)

    assert len(p.a) == len(p.b) == 2

    assert pytest.approx(p.a[0], 1e-15) == 1
    assert pytest.approx(p.a[1], 1e-15) == 7/8
    assert pytest.approx(p.b[0], 1e-15) == 1
    assert pytest.approx(p.b[1], 1e-15) == 13/8
