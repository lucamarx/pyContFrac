"""
Padé Approximant Test
"""
import math, pytest

from contfrac.pade import PadeApprox

# see https://en.wikipedia.org/wiki/Padé_table
EXP_PQ = {
    # (l,m): (a,b)
    (0,0): ([1],                      [1]),
    (0,1): ([1],                      [1,-1]),
    (0,2): ([1],                      [1,-1,1/2]),
    (0,3): ([1],                      [1,-1,1/2,-1/6]),
    (0,4): ([1],                      [1,-1,1/2,-1/6,1/24]),

    (1,0): ([1,1],                    [1]),
    (1,1): ([1,1/2],                  [1,-1/2]),
    (1,2): ([1,1/3],                  [1,-2/3,1/6]),
    (1,3): ([1,1/4],                  [1,-3/4,1/4,-1/24]),
    (1,4): ([1,1/5],                  [1,-4/5,3/10,-1/15,1/120]),

    (2,0): ([1,1,1/2],                [1]),
    (2,1): ([1,2/3,1/6],              [1,-1/3]),
    (2,2): ([1,1/2,1/12],             [1,-1/2,1/12]),
    (2,3): ([1,2/5,1/20],             [1,-3/5,3/20,-1/60]),
    (2,4): ([1,1/3,1/30],             [1,-2/3,1/5,-1/30,1/360]),

    (3,0): ([1,1,1/2,1/6],            [1]),
    (3,1): ([1,3/4,1/4,1/24],         [1,-1/4]),
    (3,2): ([1,3/5,3/20,1/60],        [1,-2/5,1/20]),
    (3,3): ([1,1/2,1/10,1/120],       [1,-1/2,1/10,-1/120]),
    (3,4): ([1,3/7,1/14,1/210],       [1,-4/7,1/7,-2/105,1/840]),

    (4,0): ([1,1,1/2,1/6,1/24],       [1]),
    (4,1): ([1,4/5,3/10,1/15,1/120],  [1,-1/5]),
    (4,2): ([1,2/3,1/5,1/30,1/360],   [1,-1/3,1/30]),
    (4,3): ([1,4/7,1/7,2/105,1/840],  [1,-3/7,1/14,-1/210]),
    (4,4): ([1,1/2,3/28,1/84,1/1680], [1,-1/2,3/28,-1/84,1/1680]),
}


def _exp_taylor():
    "exponential Taylor coefficients"
    n = 0
    while True:
        yield 1 / math.factorial(n)
        n += 1


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
    "test exponential approximant"
    for l in range(5):
        for m in range(5):
            p = PadeApprox(l, m, _exp_taylor())
            print(p)

            t = EXP_PQ[(l,m)]

            assert len(p.a) == len(t[0])

            for a,ta in zip(p.a, t[0]):
                assert pytest.approx(a, 1e-15) == ta

            assert len(p.b) == len(t[1])

            for b,tb in zip(p.b, t[1]):
                assert pytest.approx(b, 1e-15) == tb
