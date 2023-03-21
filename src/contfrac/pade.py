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

    for i in range(m):
        B[i] = -C[l+1+i]

        for j in range(m):
            k = l-m+1+i+j
            A[i,j] = C[k] if k >= 0 else 0

    # solve linear system for Q's coefficients
    b = np.flip(np.linalg.solve(A, B))

    # insert b_0=1 value
    b = np.insert(b, 0, [1])

    # compute P's coefficients
    a = np.zeros(l+1, dtype=np.float64)

    for i in range(0,l+1):
        for j in range(min(m,i)+1):
            a[i] += b[j]*C[i-j]

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


    def P(self, x : float) -> float:
        """Evaluate approximant numerator

        """
        # 1, x, x^2, ..., x^l
        x_powers = [x**i for i in range(self.l + 1)]

        # a_0 + a_1 x + a_2 x^2 + ... + a_l x^l
        return sum([a*x for a,x in zip(self.a, x_powers)])


    def Q(self, x : float) -> float:
        """Evaluate approximant denominator

        """
        # 1, x, x^2, ..., x^m
        x_powers = [x**i for i in range(self.m + 1)]

        # 1 + b_1 x + b_2 x^2 + ... + b_m x^m
        return sum([b*x for b,x in zip(self.b, x_powers)])


    def __call__(self, x : float) -> float:
        """Evaluate approximant

        """
        # 1, x, x^2, ..., x^(m+l+1)
        x_powers = [x**i for i in range(max(self.l,self.m) + 1)]

        # a_0 + a_1 x + a_2 x^2 + ... + a_l x^l
        P = sum([a*x for a,x in zip(self.a, x_powers)])

        # 1 + b_1 x + b_2 x^2 + ... + b_m x^m
        Q = sum([b*x for b,x in zip(self.b, x_powers)])

        if Q == 0: raise OverflowError

        return P / Q
