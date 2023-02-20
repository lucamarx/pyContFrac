"""
Continued fractions
"""
from __future__ import annotations

import math, numbers, itertools

from typing import Optional, Union, Callable, Generator, List


def _coefficients_from_value(x : float) -> Generator:
    while True:
        a = math.floor(x)

        yield a

        if abs(x-a) < 1e-12:
            break

        x = 1 / (x-a)


def _convergents(x : Generator) -> Generator:
    p_km1, q_km1 = 1, 0
    p_km2, q_km2 = 0, 1

    while True:
        try:
            a_k = next(x)
        except StopIteration:
            break

        p_k, q_k = a_k*p_km1 + p_km2, a_k*q_km1 + q_km2

        yield (p_k, q_k, p_k/q_k)

        p_km2, q_km2 = p_km1, q_km1
        p_km1, q_km1 = p_k, q_k


def _mobius_transform(a : int, b : int, c : int, d : int, x : Generator) -> Generator:
    while True:
        if c != 0 and d != 0 and math.floor(a/c) == math.floor(b/d):
            # emit and EGEST

            q = math.floor(a/c)

            yield q

            a, b, c, d = c, d, a-c*q, b-d*q

        else:
            # get one more term from x and INGEST
            try:
                p = next(x)
            except StopIteration:
                break

            a, b, c, d = b, a+b*p, d, c+d*p

    yield b


class ContFrac():
    """A continued fraction

    """

    def __init__(self, x : Union[numbers.Number, Callable[[], Generator]]):
        """Initialize a CF

        Parameters
        ----------

        x : Real, Callable
        A specific value or a function that returns a stream
        of coefficients

        """
        if isinstance(x, Callable):
            self._coefficients = x
        elif isinstance(x, numbers.Real):
            self._coefficients = lambda: _coefficients_from_value(float(x))
        else:
            raise TypeError


    def __str__(self) -> str:
        p, q, v = self.convergents_list()[-1]
        l = max(len(str(p)), len(str(q)))

        return f"""
  {p}
  {"-"*l} = {v}
  {q}
        """


    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.__float__()})"


    def coefficients(self) -> Generator:
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


    def coefficients_list(self, N : int = 40) -> List:
        """The coefficients as a list

        """
        return list(itertools.islice(self.coefficients(), N))


    def convergents(self) -> Generator:
        """The convergents stream

        Returns
        -------

        `[(p0, q0, p0/q0), (p1, q1, p1/q1), ...]`

        where `pn / qn` are the rational approximations to the continued
        fraction

        """
        return _convergents(self._coefficients())


    def convergents_list(self, N : int = 40) -> List:
        """The convergents as a list

        """
        return list(itertools.islice(self.convergents(), N))


    def mobius(self, a : int, b : int, c : int, d : int) -> ContFrac:
        """Apply a Möbius transform

        `       a + b·x
        ` x' = ---------
        `       c + d·x

        """
        return ContFrac(lambda: _mobius_transform(a, b, c, d, self._coefficients()))


    # UNARY ARITHMETIC

    def __int__(self) -> int:
        return int(self.__float__())


    def __float__(self) -> float:
        """Compute the best approximation up to machine precision

        """
        # TODO: improve!
        return self.convergents_list()[-1][2]


    def __trunc__(self) -> float:
        return math.trunc(self.__float__())


    def __ceil__(self) -> int:
        return math.ceil(self.__float__())


    def __floor__(self) -> int:
        return math.floor(self.__float__())


    def __round__(self, ndigits : Optional[int] = None) -> float:
        return round(self.__float__(), ndigits) if ndigits else self.__float__()


    def __pos__(self) -> ContFrac:
        return self


    def __neg__(self) -> ContFrac:
        return self.mobius(0, -1, 1, 0)


    # BINARY ARITHMETIC

    def __add__(self, other : numbers.Rational) -> ContFrac:
        n, d = other.numerator, other.denominator
        return self.mobius(n, d, d, 0)


    def __sub__(self, other : numbers.Rational) -> ContFrac:
        n, d = other.numerator, other.denominator
        return self.mobius(-n, d, d, 0)


    def __mul__(self, other : numbers.Rational) -> ContFrac:
        n, d = other.numerator, other.denominator
        return self.mobius(0, n, d, 0)


    def __truediv__(self, other : numbers.Rational) -> ContFrac:
        n, d = other.numerator, other.denominator
        return self.mobius(0, d, n, 0)
