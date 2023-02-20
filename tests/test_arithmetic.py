"""
Rational Arithmetic Test
"""
import pytest, random

from fractions import Fraction
from contfrac import ContFrac


N_ITERS = 1_000

def test_creation_from_rational():
    "create a CF from a positive rational number"
    for _ in range(N_ITERS):
        a_rt = Fraction(random.randint(1, 10_000), random.randint(1, 10_000))
        a_cf = ContFrac(a_rt)

        cv = a_cf.convergents_list()[-1]

        assert cv[0] == a_rt.numerator
        assert cv[1] == a_rt.denominator


def test_creation_from_real():
    "create a CF from a positive real number"
    for _ in range(N_ITERS):
        x = 100_000.0 * random.random()

        z = ContFrac(x)

        _, _, v = z.convergents_list()[-1]

        assert pytest.approx(v, 1e-15) == x
