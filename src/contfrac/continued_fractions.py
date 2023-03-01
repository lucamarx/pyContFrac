"""
Continued fractions
"""
from __future__ import annotations

import math, itertools, fractions

from typing import Optional, Union, Callable, Generator, List, Tuple


def _residue(x : Generator[int, None, None]) -> Generator[int, None, None]:
    """Skip first term

    """
    next(x)

    return x


def _euclid(x : fractions.Fraction) -> Generator[int, None, None]:
    """Euclid's algorithm

    """
    n, d = x.numerator, x.denominator

    while True:
        q, r = n//d, n%d

        yield q

        if r == 0: break

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


def _homographic_transform(x : Generator[int, None, None], a : int, b : int, c : int, d : int) -> Generator[int, None, None]:
    if c == 0 and d == 0:
        # ∞
        return

    while True:

        if c != 0 and d != 0 and c*d > 0 and math.trunc(a/c) == math.trunc(b/d) != 0:
            # emit next coefficient and EGEST

            q = math.trunc(a/c)

            yield q

            a, b, c, d = c, d, a-q*c, b-q*d

        else:
            # get one more term from x and INGEST
            try:
                p = next(x)

                a, b, c, d = p*a+b, a, p*c+d, c

            except StopIteration:
                if c != 0:
                    for q in _euclid(fractions.Fraction(a, c)):
                        yield q

                break


def _bihomographic_transform(x : Generator[int, None, None],
                             y : Generator[int, None, None],
                             a : int, b : int, c : int, d : int,
                             e : int, f : int, g : int, h : int) -> Generator[int, None, None]:
    if e == 0 and f == 0 and g == 0 and h == 0:
        # ∞
        return

    if e < 0 or f < 0 or g < 0 or h < 0:
        raise Exception("e, f, g, h must be positive")

    stop_x = stop_y = False

    while True:

        if stop_x and stop_y:
            # expand last term with Euclid
            if e != 0:
                for q in _euclid(fractions.Fraction(a, e)):
                    yield q
            break

        if e != 0 and f != 0 and g != 0 and h != 0 and math.trunc(a/e) == math.trunc(b/f) == math.trunc(c/g) == math.trunc(d/h) != 0:
            # emit next coefficient and EGEST

            r = math.trunc(a/e)

            yield r

            a, b, c, d, e, f, g, h = e, f, g, h, a-e*r, b-f*r, c-g*r, d-h*r

        else:
            # we need more info ...

            if (not stop_x and not stop_y and f != 0 and g != 0 and h != 0 and abs(b/f - d/h) > abs(c/g - d/h)) or stop_y:
                # get one more term from x and INGEST
                try:
                    p = next(x)

                    a, b, c, d, e, f, g, h = a*p+c, b*p+d, a, b, e*p+g, f*p+h, e, f

                except StopIteration:
                    stop_x = True

            elif not stop_y:
                # get one more term from y and INGEST
                try:
                    q = next(y)

                    a, b, c, d, e, f, g, h = a*q+b, a, c*q+d, c, e*q+f, e, g*q+h, g

                except StopIteration:
                    stop_y = True


class ContFrac():
    """A continued fraction

    """

    def __init__(self, x : Union[int, float, List[int], fractions.Fraction, Callable[[], Generator[int, None, None]]]):
        """Initialize a continued fraction

        Parameters
        ----------

        x : int, float, list, fraction, callable
        A specific value, a list or a function
        that returns a stream of coefficients

        """
        if isinstance(x, Callable):
            self._coefficients = x

        elif isinstance(x, int):
            self._coefficients = lambda: _euclid(fractions.Fraction(x,1))

        elif isinstance(x, float):
            self._coefficients = lambda: _euclid(fractions.Fraction(x))

        elif isinstance(x, List):
            self._coefficients = lambda: (a for a in x)

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


    def bihomographic(self,
                      other : ContFrac,
                      a : int, b : int, c : int, d : int,
                      e : int, f : int, g : int, h : int) -> ContFrac:
        """Apply the bihomographic transform

        `      a·x·y + b·x + c·y + d
        ` z = -----------------------
        `      e·x·y + f·x + g·y + h

        """
        return ContFrac(lambda: _bihomographic_transform(self._coefficients(),
                                                         other._coefficients(),
                                                         a, b, c, d,
                                                         e, f, g, h))


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


    def split(self) -> Tuple[int, ContFrac]:
        """Split first term from the rest

        """
        coeff = self._coefficients()

        try:
            return (next(coeff), ContFrac(lambda: _residue(self._coefficients())))
        except StopIteration:
            raise OverflowError


    # UNARY ARITHMETIC

    def __int__(self) -> int:
        try:
            coeff = self._coefficients()

            a = next(coeff)

            # CF is positive return first coefficient
            if a >= 0:
                return a

            try:
                next(coeff)

                # CF is negative
                return a+1

            except StopIteration:
                # CF is a negative integer
                return a

        except StopIteration:
            raise OverflowError


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
        return round(float(self), ndigits)


    def __pos__(self) -> ContFrac:
        return self


    def __neg__(self) -> ContFrac:
        a, r = self.split()

        return r.homographic(-a, -1, 1, 0)


    def __abs__(self) -> ContFrac:
        a, r = self.split()

        return r.homographic(-a, -1, 1, 0) if a < 0 else self


    # COMPARISON

    def __eq__(self, other : Union[int, fractions.Fraction, ContFrac]) -> bool:
        if isinstance(other, int) or isinstance(other, fractions.Fraction):
            return self.as_rational() == other

        if isinstance(other, ContFrac):
            coeff_a, coeff_b = self._coefficients(), other._coefficients()

            for a, b in zip(coeff_a, coeff_b):
                if a != b: return False

            return next(coeff_a, None) == next(coeff_b, None) == None

        return NotImplemented


    def __lt__(self, other : Union[int, fractions.Fraction, ContFrac]) -> bool:
        if isinstance(other, int) or isinstance(other, fractions.Fraction):
            return self.as_rational() < other

        if isinstance(other, ContFrac):
            coeff_a, coeff_b = self._coefficients(), other._coefficients()

            n = 0
            for a, b in zip(coeff_a, coeff_b):
                if a != b:
                    return a<b if n % 2 == 0 else b<a
                n += 1

            rest_a, rest_b = next(coeff_a, None), next(coeff_b, None)

            return (rest_a if n % 2 == 0 else rest_b) is not None

        return NotImplemented


    def __le__(self, other : Union[int, fractions.Fraction, ContFrac]) -> bool:
        return self < other or self == other


    def __gt__(self, other : Union[int, fractions.Fraction, ContFrac]) -> bool:
        return not (self < other)


    def __ge__(self, other : Union[int, fractions.Fraction, ContFrac]) -> bool:
        return self > other or self == other


    # BINARY ARITHMETIC

    def __add__(self, other : Union[int, fractions.Fraction, ContFrac]) -> ContFrac:
        if isinstance(other, int):
            return self.homographic(1, other, 0, 1)

        if isinstance(other, fractions.Fraction):
            return self.homographic(other.denominator, other.numerator, 0, other.denominator)

        if isinstance(other, ContFrac):
            return self.bihomographic(other, 0, 1, 1, 0, 0, 0, 0, 1)

        return NotImplemented


    def __radd__(self, other : Union[int, fractions.Fraction]) -> ContFrac:
        if isinstance(other, int):
            return self.homographic(1, other, 0, 1)

        if isinstance(other, fractions.Fraction):
            return self.homographic(other.denominator, other.numerator, 0, other.denominator)

        return NotImplemented


    def __sub__(self, other : Union[int, fractions.Fraction, ContFrac]) -> ContFrac:
        if isinstance(other, int):
            return self.homographic(1, -other, 0, 1)

        if isinstance(other, fractions.Fraction):
            return self.homographic(other.denominator, -other.numerator, 0, other.denominator)

        if isinstance(other, ContFrac):
            return self.bihomographic(other, 0, 1, -1, 0, 0, 0, 0, 1)

        return NotImplemented


    def __rsub__(self, other : Union[int, fractions.Fraction]) -> ContFrac:
        if isinstance(other, int):
            return self.homographic(-1, other, 0, 1)

        if isinstance(other, fractions.Fraction):
            return self.homographic(-other.denominator, other.numerator, 0, other.denominator)

        return NotImplemented


    def __mul__(self, other : Union[int, fractions.Fraction, ContFrac]) -> ContFrac:
        if isinstance(other, int):
            return self.homographic(other, 0, 0, 1)

        if isinstance(other, fractions.Fraction):
            return self.homographic(other.numerator, 0, 0, other.denominator)

        if isinstance(other, ContFrac):
            return self.bihomographic(other, 1, 0, 0, 0, 0, 0, 0, 1)

        return NotImplemented


    def __rmul__(self, other : Union[int, fractions.Fraction]) -> ContFrac:
        if isinstance(other, int):
            return self.homographic(other, 0, 0, 1)

        if isinstance(other, fractions.Fraction):
            return self.homographic(other.numerator, 0, 0, other.denominator)

        return NotImplemented


    def __truediv__(self, other : Union[int, fractions.Fraction, ContFrac]) -> ContFrac:
        if isinstance(other, int) or isinstance(other, fractions.Fraction):

            if isinstance(other, int):
                n, d = other, 1
            else:
                n, d = other.numerator, other.denominator

                if d < 0: n, d = -n, -d

            if n == 0: raise OverflowError

            return self.homographic(d, 0, 0, n)

        if isinstance(other, ContFrac):
            if other == 0: raise OverflowError

            if other > 0:
                return self.bihomographic(other, 0, 1, 0, 0, 0, 0, 1, 0)
            else:
                return self.bihomographic(-other, 0, -1, 0, 0, 0, 0, 1, 0)


    def __rtruediv__(self, other : Union[int, fractions.Fraction]) -> ContFrac:
        if isinstance(other, int):
            n, d = other, 1
        elif isinstance(other, fractions.Fraction):
            n, d = other.numerator, other.denominator

            if d < 0: n, d = -n, -d
        else:
            return NotImplemented

        if self == 0: raise OverflowError

        return self.homographic(0, n, d, 0)
