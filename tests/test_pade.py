"""
Padé Approximant Test
"""
import pytest

from contfrac.pade import PadeApprox

EXP_PQ = {
    # (l,m): (a,b)
    (0,0): ([1],          [1]),
    (0,1): ([1],          [1,-1]),
    (0,2): ([1],          [1,-1,1/2]),

    (1,0): ([1,1],        [1]),
    (1,1): ([1,1/2],      [1,-1/2]),
    (1,2): ([1,1/3],      [1,-2/3,1/6]),

    (2,0): ([1,1,1/2],    [1]),
    (2,1): ([1,2/3,1/6],  [1,-1/3]),
    (2,2): ([1,1/2,1/12], [1,-1/2,1/12])
}

def test_basic():
    "test basic function approximation"
    c = [1,-3/4,39/32]

    p = PadeApprox(1, 1, c)

    assert len(p.a) == len(p.b) == 2

    assert pytest.approx(p.a[0], 1e-15) == 1
    assert pytest.approx(p.a[1], 1e-15) == 7/8
    assert pytest.approx(p.b[0], 1e-15) == 1
    assert pytest.approx(p.b[1], 1e-15) == 13/8


def test_exponential():
    "test exponential approximation"
    c = [1, 1, 1/2, 1/6, 1/24, 1/120]

    for l in range(3):
        for m in range(3):
            p = PadeApprox(l, m, c)
            print(p)

            t = EXP_PQ[(l,m)]

            assert len(p.a) == len(t[0])

            for a,ta in zip(p.a, t[0]):
                assert pytest.approx(a, 1e-15) == ta

            assert len(p.b) == len(t[1])

            for b,tb in zip(p.b, t[1]):
                assert pytest.approx(b, 1e-15) == tb
