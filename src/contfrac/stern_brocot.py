"""
Stern-Brocot numerals
"""
from __future__ import annotations

import numpy as np
import re, fractions

from typing import Union


L = np.array([[1, 1],
              [0, 1]],
             dtype=np.int64)


R = np.array([[1, 0],
              [1, 1]],
             dtype=np.int64)


def _mediant(S : np.ndarray) -> fractions.Fraction:
    """Compute the mediant

    """
    return fractions.Fraction(S[1,0] + S[1,1], S[0,0] + S[0,1])


def decode(S : str) -> fractions.Fraction:
    """Decode a Stern-Brocot numeral into a fraction

    """
    T = np.array([[1, 0],
                  [0, 1]],
                 dtype=np.int64)

    for s in S:
        if s == "L":
            T = T @ L
        elif s == "R":
            T = T @ R
        else:
            raise ValueError

    return _mediant(T)


def encode(r : fractions.Fraction) -> str:
    """Encode a fraction into a Stern-Brocot numeral

    """
    S, T = "", np.array([[1, 0],
                         [0, 1]],
                        dtype=np.int64)

    while r != _mediant(T):
        if r < _mediant(T):
            S, T = S + "L", T @ L
        else:
            S, T = S + "R", T @ R

    return S


def homographic(S : str, a :int, b : int, c : int, d : int) -> str:
    """Apply the homographic transform to the Stern-Brocot numeral

    """
    T = np.array([[c, d],
                  [a, b]],
                 dtype=np.int64)

    for s in S:
        if s == "L":
            T = T @ L
        elif s == "R":
            T = T @ R
        else:
            raise ValueError

    return encode(_mediant(T))


class SterBroc():
    """Signed Stern-Brocot numbers

    """

    def __init__(self, x : Union[str, fractions.Fraction]):
        """Initialize a Stern-Brocot number

        """
        if isinstance(x, str):
            m = re.match("^([+-]*)([LR]*)$", x)
            if m:
                self.sign = -1 if m.group(1) == "-" else +1
                self.value = decode(m.group(2))
                self.numeral = m.group(2)
            else:
                raise Exception("invalid Stern-Brocot numeral")

        elif isinstance(x, fractions.Fraction):
            if x == 0:
                self.sign = 0
                self.value = 0
                self.numeral = None
            else:
                self.sign = -1 if x < 0 else 1
                self.value = abs(x)
                self.numeral = encode(abs(x))

        else:
            raise TypeError


    def __str__(self) -> str:
        if self.sign == 0:
            return "0"

        s = "+" if self.sign == 1 else "-"
        n, d = str(self.value.numerator), str(self.value.denominator)
        l = max(len(n), len(d))

        return f"""
    {n}
  {s} {'-'*l} = {s}{self.numeral}
    {d}
        """


    def __repr__(self) -> str:
        if self.sign == 0:
            return f"{self.__class__.__name__}(0)"

        s = "+" if self.sign == 1 else "-"

        return f"{self.__class__.__name__}({s}{self.numeral})"


    def homographic(self, a : int, b : int, c : int, d : int) -> SterBroc:
        """Apply the homographic transform

        `      a·x + b
        ` y = ---------
        `      c·x + d

        """
        pass
