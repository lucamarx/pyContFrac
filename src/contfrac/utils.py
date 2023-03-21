"""Utilities

"""
import re
import numpy as np
from typing import Generator


def latex_format_number(x : float) -> str:
    """Format a number in latex

    """
    m = re.match("^(.+)E(.+)$", f"{x:.2E}")

    if m is None:
        return str(x)

    if m.group(2) == "+00":
        return f"{m.group(1)}"

    return f"{m.group(1)} \\cdot 10^{{{m.group(2)}}}"


def latex_format_poly(a : np.ndarray, max_order : int = 3) -> str:
    """Format polynomial in latex

    """
    powers = [f"z^{i}" if i>0 else "" for i in range(len(a))]
    terms = [f"{latex_format_number(a[i])} {powers[i]}" for i in range(min(len(a), max_order))]

    if len(terms) < len(a):
       terms.append("\\cdots")
       terms.append(f"{latex_format_number(a[-1])} {powers[-1]}")

    return " + ".join(terms)


class CachedGenerator():
    """A generator wrapper that remembers the last values

    """

    def __init__(self, iter : Generator, size : int = 10):
        self.iter = iter
        self.size = size
        self.last_values = []


    def __repr__(self):
        return f"{self.__class__.__name__}({self.size})"


    def __iter__(self):
        return self


    def __next__(self):
        v = next(self.iter, None)

        self.last_values = [v] + self.last_values

        if len(self.last_values) > self.size:
            self.last_values.pop()

        if v is None: raise StopIteration

        return v
