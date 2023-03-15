"""
Stern-Brocot numerals
"""
import numpy as np
import fractions


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
