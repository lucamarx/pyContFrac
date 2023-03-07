"""
Utilities
"""
from typing import Generator


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
