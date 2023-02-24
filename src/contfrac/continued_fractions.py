"""
Continued fractions
"""
from __future__ import annotations

import math, itertools, fractions

from typing import Optional, Union, Callable, Generator, List, Tuple


def _euclid(x : fractions.Fraction) -> Generator[int, None, None]:
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



def _ingest(p : int, a : int, b : int, c : int, d : int) -> Tuple[int, int, int, int]:
    return (p*a+b, a, p*c+d, c)


def _egest(q : int, a : int, b : int, c : int, d : int) -> Tuple[int, int, int, int]:
    return (c, d, a-q*c, b-q*d)


def _ingest_egest_check(a : int, b : int, c : int, d : int) -> Optional[int]:
    if c == 0 or d == 0:
        return None

    f1, f2 = math.floor(a/c), math.floor(b/d)

    if f1 == f2 or abs(f1-f2) == 1:
        q = min(f1, f2)

        if q != 0: return q

    return None


def _homographic_transform(x : Generator[int, None, None], a : int, b : int, c : int, d : int) -> Generator[int, None, None]:
    while True:
        if c == 0:
            break

        q = _ingest_egest_check(a, b, c, d)

        if q is not None:
            # emit next coefficient and EGEST

            yield q

            a, b, c, d = _egest(q, a, b, c, d)

        else:
            # get one more term from x and INGEST
            try:
                p = next(x)

                a, b, c, d = _ingest(p, a, b, c, d)

            except StopIteration:
                # TODO: go on expanding
                break


class ContFrac():
    """A continued fraction

    """

    def __init__(self, x : Union[float, fractions.Fraction, Callable[[], Generator[int, None, None]]]):
        """Initialize a CF

        Parameters
        ----------

        x : Real, Rational, Callable
        A specific value or a function that returns a stream
        of coefficients

        """
        if isinstance(x, Callable):
            self._coefficients = x
        elif isinstance(x, float):
            self._coefficients = lambda: _euclid(fractions.Fraction(x))
        elif isinstance(x, fractions.Fraction):
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

        `      a·x + b
        ` z = ---------
        `      c·x + d

        """
        return ContFrac(lambda: _homographic_transform(self._coefficients(), a, b, c, d))


    def as_rational(self) -> fractions.Fraction:
        """Compute the best rational approximation up to machine precision, it
        is exact if the number was rational in the first place

        """
        best_conv = None
        for conv in self.convergents():
            if best_conv is not None and abs(conv - best_conv) < 1e-14:
                break

            best_conv = conv

        # ∞
        if best_conv is None: raise OverflowError

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

    def __add__(self, other : fractions.Fraction) -> ContFrac:
        n, d = other.numerator, other.denominator
        return self.homographic(d, n, 0, d)


    def __sub__(self, other : fractions.Fraction) -> ContFrac:
        n, d = other.numerator, other.denominator
        return self.homographic(d, -n, 0, d)


    def __mul__(self, other : fractions.Fraction) -> ContFrac:
        n, d = other.numerator, other.denominator
        return self.homographic(n, 0, 0, d)


    def __truediv__(self, other : fractions.Fraction) -> ContFrac:
        n, d = other.numerator, other.denominator

        if n == 0: raise OverflowError

        return self.homographic(d, 0, 0, n)
