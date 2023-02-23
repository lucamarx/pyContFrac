"""
Continued fractions
"""
from __future__ import annotations

import math, numbers, itertools, fractions

from typing import Optional, Union, Callable, Generator, List, Tuple


def _euclid(x : numbers.Rational) -> Generator[int, None, None]:
    """Euclid's algorithm

    """
    n, d = x.numerator, x.denominator
    while True:
        q, r = n//d, n%d

        yield q

        if r == 0:
            break

        n, d = d, r


def _convergents(x : Generator) -> Generator[fractions.Fraction, None, None]:
    p_km1, q_km1 = 1, 0
    p_km2, q_km2 = 0, 1

    while True:
        try:
            a_k = next(x)
        except StopIteration:
            break

        p_k, q_k = a_k*p_km1 + p_km2, a_k*q_km1 + q_km2

        yield fractions.Fraction(p_k, q_k)

        p_km2, q_km2 = p_km1, q_km1
        p_km1, q_km1 = p_k, q_k


def _homographic_transform(a : int, b : int, c : int, d : int, x : Generator[int, None, None]) -> Generator[int, None, None]:
    while True:
        if c == 0 and d == 0:
            break

        elif c != 0 and d != 0 and math.floor(a/c) == math.floor(b/d):
            # emit next coefficient and EGEST

            q = math.floor(a/c)

            yield q

            a, b, c, d = c, d, a-c*q, b-d*q

        else:
            # get one more term from x and INGEST
            try:
                p = next(x)

                # p ≠ ∞
                a, b, c, d = b, a+b*p, d, c+d*p

            except StopIteration:
                # p = ∞
                a, b, c, d = b, b, d, d


class ContFrac():
    """A continued fraction

    """

    def __init__(self, x : Union[numbers.Real, numbers.Rational, Callable[[], Generator[int, None, None]]]):
        """Initialize a CF

        Parameters
        ----------

        x : Real, Rational, Callable
        A specific value or a function that returns a stream
        of coefficients

        """
        if isinstance(x, Callable):
            self._coefficients = x
        elif isinstance(x, numbers.Real):
            self._coefficients = lambda: _euclid(fractions.Fraction(*x.as_integer_ratio()))
        elif isinstance(x, numbers.Rational):
            self._coefficients = lambda: _euclid(x)
        else:
            raise TypeError


    def __str__(self) -> str:
        v = self.as_rational()
        l = max(len(str(v.numerator)), len(str(v.denominator)))

        return f"""
  {v.numerator}
  {"-" * l} = {v}
  {v.denominator}
        """


    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.as_rational()})"


    def coefficients(self) -> Generator[int, None, None]:
        """The coefficients stream of the continued fraction

        `            1
        ` a0 + -------------
        `              1
        `      a1 + --------
        `                 1
        `           a2 + ---
        `                ...

        Returns
        -------

        `[a0, a1, a2, ...]`

        """
        return self._coefficients()


    def coefficients_as_list(self, N : int = 40) -> List[int]:
        """The coefficients as a list

        """
        return list(itertools.islice(self.coefficients(), N))


    def convergents(self) -> Generator[fractions.Fraction, None, None]:
        """The convergents stream

        Returns
        -------

        `[p0/q0, p1/q1, ...]`

        where `pn / qn` are the rational approximations to the continued
        fraction

        """
        return _convergents(self._coefficients())


    def convergents_as_list(self, N : int = 40) -> List[fractions.Fraction]:
        """The convergents as a list

        """
        return list(itertools.islice(self.convergents(), N))


    def homographic(self, a : int, b : int, c : int, d : int) -> ContFrac:
        """Apply the homographic transform

        `      a + b·x
        ` z = ---------
        `      c + d·x

        """
        return ContFrac(lambda: _homographic_transform(a, b, c, d, self._coefficients()))


    def as_rational(self) -> fractions.Fraction:
        """Compute the best rational approximation up to machine precision, it
        is exact if the number was rational in the first place

        """
        best_conv = None
        for conv in self.convergents():
            if best_conv is not None and abs(conv - best_conv) < 1e-14:
                break

            best_conv = conv

        if best_conv is None:
            # ∞
            raise OverflowError

        return best_conv


    # UNARY ARITHMETIC

    def __int__(self) -> int:
        return int(self.as_rational())


    def __float__(self) -> float:
        return float(self.as_rational())


    def as_integer_ratio(self) -> Tuple[int, int]:
        return self.as_rational().as_integer_ratio()


    def __trunc__(self) -> float:
        return math.trunc(self.as_rational())


    def __ceil__(self) -> int:
        return math.ceil(self.as_rational())


    def __floor__(self) -> int:
        return math.floor(self.as_rational())


    def __round__(self, ndigits : Optional[int] = None) -> float:
        return round(float(self), ndigits) if ndigits else float(self)


    def __pos__(self) -> ContFrac:
        return self


    def __neg__(self) -> ContFrac:
        return self.homographic(0, -1, 1, 0)


    # BINARY ARITHMETIC

    def __add__(self, other : numbers.Rational) -> ContFrac:
        n, d = other.numerator, other.denominator
        return self.homographic(n, d, d, 0)


    def __sub__(self, other : numbers.Rational) -> ContFrac:
        n, d = other.numerator, other.denominator
        return self.homographic(-n, d, d, 0)


    def __mul__(self, other : numbers.Rational) -> ContFrac:
        n, d = other.numerator, other.denominator
        return self.homographic(0, n, d, 0)


    def __truediv__(self, other : numbers.Rational) -> ContFrac:
        n, d = other.numerator, other.denominator
        return self.homographic(0, d, n, 0)
