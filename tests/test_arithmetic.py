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
        r = Fraction(random.randint(1, 10_000), random.randint(1, 10_000))

        assert ContFrac(r).as_rational() == r


def test_creation_from_real():
    "create a CF from a positive real number"
    for _ in range(N_ITERS):
        x = 100_000.0 * random.random()

        assert pytest.approx(float(ContFrac(x)), 1e-15) == x


def test_homographic():
    "test homographic transform"
    for _ in range(N_ITERS):
        a, b, c, d = random.randint(-100, 100), random.randint(-100, 100), random.randint(1, 100), random.randint(1, 100)

        r = Fraction(random.randint(1, 100), random.randint(1, 100))

        s = ContFrac(r).homographic(a, b, c, d).as_rational()

        assert ((a*r + b) / (c*r + d)) == s


