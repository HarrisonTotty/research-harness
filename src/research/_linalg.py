"""Exact linear algebra over the rationals and prime fields.

Representability questions are only meaningful relative to a field (Matroid
page, terminology hazard), so the elimination here is generic over an
:class:`ExactField` — a trio of arithmetic hooks — with instances for the
rationals and for ``GF(p)``. :mod:`research.matroid` uses the column-rank
routine to build linear matroids and :mod:`research.positroid` the exact
determinant to test the sign of maximal minors.
"""

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from fractions import Fraction

__all__ = [
    "RATIONALS",
    "ExactField",
    "column_rank",
    "det_q",
    "gf_scalar",
    "is_prime",
    "linear_independence_family",
    "prime_field",
]

_SMALLEST_PRIME: int = 2
"""Smallest field characteristic accepted by :func:`prime_field`."""


def is_prime(n: int) -> bool:
    """Return whether ``n`` is a prime number (trial division; small inputs)."""
    if n < _SMALLEST_PRIME:
        return False
    return all(n % d != 0 for d in range(2, math.isqrt(n) + 1))


def gf_scalar(value: Fraction | int, p: int) -> int:
    """Reduce an exact scalar into ``GF(p)``.

    A fraction ``a/b`` maps to ``a * b^(p-2) mod p`` — Fermat inversion of
    the denominator, valid because ``p`` is prime.

    Raises:
        ValueError: If the denominator is divisible by ``p``, so the scalar
            has no image in ``GF(p)``.
    """
    frac = Fraction(value)
    if frac.denominator % p == 0:
        msg = (
            f"vector entry {value} has denominator divisible by {p}, "
            f"so it has no image in GF({p})"
        )
        raise ValueError(msg)
    return frac.numerator * pow(frac.denominator, p - 2, p) % p


@dataclass(frozen=True)
class ExactField[S]:
    """Arithmetic hooks driving exact Gaussian elimination over one field.

    Elements must be canonical — one representative per field element, with
    the zero element falsy — so that elimination can test for zero by truth
    value. ``sub_mul(a, c, b)`` returns ``a - c * b``, ``inv`` the inverse
    of a nonzero element, and ``mul`` the product.
    """

    sub_mul: Callable[[S, S, S], S]
    inv: Callable[[S], S]
    mul: Callable[[S, S], S]


RATIONALS: ExactField[Fraction] = ExactField(
    sub_mul=lambda a, c, b: a - c * b,
    inv=lambda a: 1 / a,
    mul=lambda a, b: a * b,
)
"""The rationals, on :class:`fractions.Fraction` representatives."""


def prime_field(p: int) -> ExactField[int]:
    """Return ``GF(p)`` arithmetic on the canonical representatives ``0..p-1``.

    Inverses use Fermat's little theorem, valid because ``p`` is prime; the
    caller is responsible for reducing inputs, e.g. via :func:`gf_scalar`.
    """
    return ExactField(
        sub_mul=lambda a, c, b: (a - c * b) % p,
        inv=lambda a: pow(a, p - 2, p),
        mul=lambda a, b: a * b % p,
    )


def column_rank[S](columns: Sequence[Sequence[S]], field: ExactField[S]) -> int:
    """Return the rank of the given column vectors over ``field``.

    Exact Gaussian elimination: each column is reduced against the pivots
    found so far and either contributes a new pivot or is dependent.
    """
    pivots: list[int] = []
    reduced: list[list[S]] = []
    for column in columns:
        vector = list(column)
        for position, pivot_vector in zip(pivots, reduced, strict=True):
            coefficient = vector[position]
            if coefficient:
                vector = [
                    field.sub_mul(a, coefficient, b)
                    for a, b in zip(vector, pivot_vector, strict=True)
                ]
        position = next((i for i, a in enumerate(vector) if a), -1)
        if position < 0:
            continue
        inverse = field.inv(vector[position])
        pivots.append(position)
        reduced.append([field.mul(a, inverse) for a in vector])
    return len(pivots)


def linear_independence_family[S](
    columns: Sequence[Sequence[S]], field: ExactField[S]
) -> frozenset[int]:
    """Return the linear independence family of the columns, as bitmasks.

    Bit ``i`` of a mask selects ``columns[i]``; a mask is in the family when
    its columns are linearly independent over ``field``. Exponential in the
    number of columns, matching the explicit design of the matroid module.
    """
    n = len(columns)
    family: set[int] = set()
    for mask in range(1 << n):
        chosen = [columns[b] for b in range(n) if mask >> b & 1]
        if column_rank(chosen, field) == mask.bit_count():
            family.add(mask)
    return frozenset(family)


def det_q(rows: Sequence[Sequence[Fraction]]) -> Fraction:
    """Return the determinant of a square matrix over the rationals.

    Exact fraction-preserving Gaussian elimination with row pivoting; the
    empty matrix has determinant one.
    """
    d = len(rows)
    matrix = [list(row) for row in rows]
    det = Fraction(1)
    for col in range(d):
        pivot_row = next((r for r in range(col, d) if matrix[r][col]), None)
        if pivot_row is None:
            return Fraction(0)
        if pivot_row != col:
            matrix[col], matrix[pivot_row] = matrix[pivot_row], matrix[col]
            det = -det
        pivot = matrix[col][col]
        det *= pivot
        inverse = 1 / pivot
        for r in range(col + 1, d):
            factor = matrix[r][col] * inverse
            if factor:
                matrix[r] = [
                    a - factor * b for a, b in zip(matrix[r], matrix[col], strict=True)
                ]
    return det
