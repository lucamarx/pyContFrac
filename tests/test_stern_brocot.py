"""
Rational Arithmetic Test
"""
import random

from fractions import Fraction
from contfrac.stern_brocot import encode, decode

N_ITERS = 10_000

STERN_BROCOT = [
    (1, 1, ""),

    (1, 2, "L"),
    (2, 1, "R"),

    (1, 3, "LL"),
    (2, 3, "LR"),
    (3, 2, "RL"),
    (3, 1, "RR"),

    (1, 4, "LLL"),
    (2, 5, "LLR"),
    (3, 5, "LRL"),
    (3, 4, "LRR"),
    (4, 3, "RLL"),
    (5, 3, "RLR"),
    (5, 2, "RRL"),
    (4, 1, "RRR"),

    (1, 5, "LLLL"),
    (2, 7, "LLLR"),
    (3, 8, "LLRL"),
    (3, 7, "LLRR"),
    (4, 7, "LRLL"),
    (5, 8, "LRLR"),
    (5, 7, "LRRL"),
    (4, 5, "LRRR"),
    (5, 4, "RLLL"),
    (7, 5, "RLLR"),
    (8, 5, "RLRL"),
    (7, 4, "RLRR"),
    (7, 3, "RRLL"),
    (8, 3, "RRLR"),
    (7, 2, "RRRL"),
    (5, 1, "RRRR"),
]


def test_encoding():
    "test encoding fractions into Stern-Brocot numerals"

    for p, q, S in STERN_BROCOT:
        assert encode(Fraction(p,q)) == S


def test_decoding():
    "test decoding Stern-Brocot numerals into fractions"

    for p, q, S in STERN_BROCOT:
        assert Fraction(p,q) == decode(S)


def test_random_decoding_encoding():
    "test decoding-encoding with random input"

    for _ in range(N_ITERS):
        S = "".join([["L", "R"][random.randint(0,1)] for _ in range(random.randint(0,10))])

        assert encode(decode(S)) == S


def test_random_encoding_decoding():
    "test encoding-decoding with random input"

    for _ in range(N_ITERS):
        r = Fraction(random.randint(1,1000), random.randint(1,1000))

        assert decode(encode(r)) == r
