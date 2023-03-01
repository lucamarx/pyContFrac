"""
Rational Arithmetic Test
"""
import math, pytest, random

from fractions import Fraction
from contfrac import ContFrac


N_ITERS = 10_000

def test_creation_from_int():
    "create a CF from a positive integer"
    for _ in range(N_ITERS):
        r = random.randint(-10_000, 10_000)

        s = ContFrac(r)

        assert s.as_rational() == Fraction(r, 1)
        assert s.as_integer_ratio() == (r, 1)

        assert int(s) == r
        assert float(s) == float(r)

        assert math.floor(s) == r
        assert math.ceil(s) == r

        assert abs(s) == abs(r)
        assert -s == -r


def test_creation_from_rational():
    "create a CF from a positive rational number"
    for _ in range(N_ITERS):
        r = Fraction(random.randint(-10_000, 10_000), random.randint(1, 10_000))

        s = ContFrac(r)

        assert s.as_rational() == r
        assert s.as_integer_ratio() == r.as_integer_ratio()

        assert int(s) == int(r)
        assert float(s) == float(r)

        assert math.floor(s) == math.floor(r)
        assert math.ceil(s) == math.ceil(r)

        assert abs(s) == abs(r)
        assert -s == -r


def test_creation_from_real():
    "create a CF from a positive real number"
    for _ in range(N_ITERS):
        x = 100_000.0 * random.random()

        y = ContFrac(x)

        assert pytest.approx(float(ContFrac(x)), 1e-15) == x

        assert int(y) == int(x)

        assert math.floor(y) == math.floor(x)
        assert math.ceil(y) == math.ceil(x)

        assert math.trunc(y) == math.trunc(x)

        for n in [None, 1, 2, 3, 4]:
            assert round(y, n) == round(x, n)


def test_creation_from_list():
    "create a CF from its coefficients"
    for _ in range(N_ITERS):
        r = [random.randint(1,10) for _ in range(random.randint(1, 10))]

        s = ContFrac(r)

        assert s.coefficients_as_list() == r


def test_comparison():
    "test comparison operators"
    for _ in range(N_ITERS):
        a = Fraction(random.randint(1, 1000), random.randint(1, 1000))
        b = Fraction(random.randint(1, 1000), random.randint(1, 1000))

        assert ContFrac(a) == a
        assert (ContFrac(a) == b) == (a == b)

        assert (ContFrac(a) <= b) == (a <= b)
        assert (ContFrac(a) < b)  == (a < b)

        assert (ContFrac(a) >= b) == (a >= b)
        assert (ContFrac(a) > b)  == (a > b)


def test_homographic():
    "test homographic transform"
    for _ in range(N_ITERS):
        a, b, c, d = random.randint(-100, 100), random.randint(-100, 100), random.randint(-100, 100), random.randint(-100, 100)

        r = Fraction(random.randint(1, 100), random.randint(1, 100))

        s = ContFrac(r).homographic(a, b, c, d)

        if (c == 0 and d == 0) or c*r + d == 0:
            # s = ∞
            assert len(s.coefficients_as_list()) == 0
        else:
            # s < ∞
            assert ((a*r + b) / (c*r + d)) == s.as_rational()


def test_bihomographic():
    "test bihomographic transform"
    for _ in range(N_ITERS):
        a, b, c, d = random.randint(-100, 100), random.randint(-100, 100), random.randint(-100, 100), random.randint(-100, 100)
        e, f, g, h = random.randint(1, 100), random.randint(1, 100), random.randint(1, 100), random.randint(1, 100)

        x = Fraction(random.randint(1, 100), random.randint(1, 100))
        y = Fraction(random.randint(1, 100), random.randint(1, 100))

        z = ContFrac(x).bihomographic(ContFrac(y), a, b, c, d, e, f, g, h)

        assert ((a*x*y + b*x + c*y + d) / (e*x*y + f*x + g*y + h)) == z.as_rational()


def test_sum():
    "test sum"
    for _ in range(N_ITERS):
        a = Fraction(random.randint(-1000, 1000), random.randint(1, 1000))
        b = Fraction(random.randint(-1000, 1000), random.randint(1, 1000))

        assert ContFrac(a) + b == a + b
        assert a + ContFrac(b) == a + b
        assert ContFrac(a) + ContFrac(b) == a + b

        c = random.randint(1, 1000)

        assert ContFrac(a) + c == a + c
        assert a + ContFrac(c) == a + c
        assert ContFrac(a) + ContFrac(c) == a + c


def test_sub():
    "test subtraction"
    for _ in range(N_ITERS):
        a = Fraction(random.randint(1, 1000), random.randint(1, 1000))
        b = Fraction(random.randint(1, 1000), random.randint(1, 1000))

        assert ContFrac(a) - b == a - b
        assert a - ContFrac(b) == a - b
        assert ContFrac(a) - ContFrac(b) == a - b

        c = random.randint(1, 1000)

        assert ContFrac(a) - c == a - c
        assert a - ContFrac(c) == a - c
        assert ContFrac(a) - ContFrac(c) == a - c


def test_mul():
    "test multiplication"
    for _ in range(N_ITERS):
        a = Fraction(random.randint(-1000, 1000), random.randint(1, 1000))
        b = Fraction(random.randint(1, 1000), random.randint(1, 1000))

        assert ContFrac(a) * b == a * b
        assert a * ContFrac(b) == a * b
        assert ContFrac(a) * ContFrac(b) == a * b

        c = random.randint(1, 1000)

        assert ContFrac(a) * c == a * c
        assert a * ContFrac(c) == a * c
        assert ContFrac(a) * ContFrac(c) == a * c


def test_div():
    "test division"
    for _ in range(N_ITERS):
        a = Fraction(random.randint(-1000, 1000), random.randint(1, 1000))
        b = Fraction(random.randint(1, 1000), random.randint(1, 1000))

        if b != 0:
            assert ContFrac(a) / b == a / b
            assert a / ContFrac(b) == a / b
            assert ContFrac(a) / ContFrac(b) == a / b

        c = random.randint(1, 1000)

        assert ContFrac(a) / c == a / c
        assert a / ContFrac(c) == a / c
        assert ContFrac(a) / ContFrac(c) == a / c
