"""
Rational Arithmetic Test
"""
import pytest, random

from fractions import Fraction
from contfrac import ContFrac


N_ITERS = 10_000

def test_creation_from_rational():
    "create a CF from a positive rational number"
    for _ in range(N_ITERS):
        a_rt = Fraction(random.randint(1, 10_000), random.randint(1, 10_000))
        a_cf = ContFrac(a_rt)

        cv = a_cf.convergents_list()[-1]

        assert cv.numerator == a_rt.numerator
        assert cv.denominator == a_rt.denominator


def test_creation_from_real():
    "create a CF from a positive real number"
    for _ in range(N_ITERS):
        x = 100_000.0 * random.random()

        z = ContFrac(x)

        cv = z.convergents_list()[-1]

        assert pytest.approx(cv, 1e-15) == x


def test_mobius():
    "test Möbius transform"
    for _ in range(N_ITERS):
        a_rt = Fraction(random.randint(1, 1000), random.randint(1, 1000))
        a_cf = ContFrac(a_rt)

        a_cv = a_cf.convergents_list()[-1]

        a, b, c, d = random.randint(1, 1000), random.randint(1, 1000), random.randint(1, 1000), random.randint(1, 1000)

        b_cf = a_cf.mobius(a, b, c, d)
        b_cv = b_cf.convergents_list()[-1]

        assert Fraction((a + b*a_cv), (c + d*a_cv)) == b_cv
