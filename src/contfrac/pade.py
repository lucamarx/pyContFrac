"""Padé Approximants

"""
import numpy as np
import itertools
from typing import Generator, List, Tuple, Union

from .utils import latex_format_poly


def _pade(l : int, m : int, c : List[float]) -> Tuple[np.ndarray, np.ndarray]:
    """Compute Padé approximant given McLaurin coefficients

    `             P(z)
    ` [l/m] (z) = ----
    `      f      Q(z)

    """
    if len(c) < l+m+1:
        raise Exception(f"f must provide at least {l+m+1} coefficients")

    C = np.array(c, dtype=np.float64)

    # setup linear system
    A = np.zeros((m, m), dtype=np.float64)
    B = np.zeros(m, dtype=np.float64)

    for j in range(m):
        B[j] = -C[l+j+1]

        for k in range(m):
            A[k,j] = C[l-m-j+k+1]

    # solve linear system for Q's coefficients
    b = np.flip(np.linalg.solve(A, B))

    # compute P's coefficients
    a, n = C[0:l+1], min(l,m)

    if n > 0:
        for i in range(1,l+1):
            a[i] += np.dot(b[0:n], np.flip(c[0:n]))

    # insert b_0=1 value
    b = np.insert(b, 0, [1])

    return (a, b)


class PadeApprox():
    """Padé approximant

    """

    def __init__(self, l : int, m : int, f : Union[List[float], Generator[float, None, None]]):
        """Initialize the l/m Padé approximant

        `             P(z)
        ` [l/m] (z) = ----
        `      f      Q(z)

        Parameters
        ----------

        l,m : int
        Approximant order

        f : generator or list
        McLaurin coefficients

        """
        self.l, self.m = l, m

        if isinstance(f, List):
            coeffs = f
        elif isinstance(f, Generator):
            coeffs = list(itertools.islice(f, l+m+1))
        else:
            raise TypeError

        self.a, self.b = _pade(l, m, coeffs[0:l+m+1])


    def __str__(self) -> str:
        return f"[{self.l}/{self.m}](z)"


    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.l}, {self.m})"


    def _repr_latex_(self) -> str:
        return f"$[{self.l}/{self.m}](x) = \\frac{{{latex_format_poly(self.a)}}}{{{latex_format_poly(self.b)}}}$"


    def __call__(self, x : float) -> float:
        """Evaluate Approximant

        """
        # 1, x, x^2, ..., x^(m+l+1)
        x_powers = [x**i for i in range(max(self.l,self.m) + 1)]

        # a_0 + a_1 x + a_2 x^2 + ... + a_l x^l
        P = sum([a*x for a,x in zip(self.a, x_powers)])

        # 1 + b_1 x + b_2 x^2 + ... + b_m x^m
        Q = sum([b*x for b,x in zip(self.b, x_powers)])

        if Q == 0: raise OverflowError

        return P / Q
