"""
Padé Approximants
"""
import numpy as np
import itertools
from typing import Generator, List, Tuple, Union


def _pade(l : int, m : int, c : List[float]) -> Tuple[np.ndarray, np.ndarray]:
    """Compute Padé approximant given McLaurin coefficients

    `             P(z)
    ` [l/m] (z) = ----
    `      f      Q(z)

    """
    C = np.array(c, dtype=np.float64)

    # setup linear system
    A = np.zeros((m, m), dtype=np.float64)
    B = np.zeros(m, dtype=np.float64)

    for j in range(m):
        B[j] = -C[l+j+1]

        for k in range(m):
            A[k,j] = C[l-m-j+k+1]

    # solve linear system for Q's coefficients
    b = np.linalg.solve(A, B)
    b = np.flip(b)

    # compute P's coefficients
    a, n = C[0:l+1], min(l,m)

    if n > 0:
        for i in range(1,l+1):
            a[i] += np.dot(b[0:n], np.flip(c[0:n]))

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
            coeffs = list(itertools.islice(f, l+m))
        else:
            raise TypeError

        if len(coeffs) < l+m+1:
            raise Exception(f"f must provide at least {l+m+1} coefficients")

        self.a, self.b = _pade(l, m, coeffs[0:(l+m+1)])


    def __str__(self) -> str:
        return f"""
  [{self.l}/{self.m}] (z)
           f
        """

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.l}, {self.m}, f)"
