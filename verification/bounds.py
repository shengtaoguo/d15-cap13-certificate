"""Exact bound formulas used by the cap-13 certificate checker."""

from __future__ import annotations


def ceil_div(numerator: int, denominator: int) -> int:
    if denominator <= 0:
        raise ValueError("nonpositive divisor")
    return -(-numerator // denominator)


def vertex_bounds(q: int, maximum_dimension: int) -> tuple[int, ...]:
    """Return the recursively propagated vertex floors in ranks 3 through d."""
    if isinstance(q, bool) or not isinstance(q, int) or q < 5:
        raise ValueError("q must be an integer at least 5")
    if (
        isinstance(maximum_dimension, bool)
        or not isinstance(maximum_dimension, int)
        or maximum_dimension < 3
    ):
        raise ValueError("maximum_dimension must be an integer at least 3")
    n_q = ceil_div(q + 4, 2)
    values = [n_q]
    for rank in range(4, maximum_dimension + 1):
        previous = values[-1]
        small_face_cap = min(5, 2 + (2 * (previous - 2)) // q)
        extension = ceil_div(q - small_face_cap + 1, 2)
        t = rank - 3
        star_bound = min(
            n_q + t * (n_q - 3),
            q + 2 + t * (n_q - 4),
            ceil_div(3 * q, 2) + 2 + t * (n_q - 5),
        )
        values.append(max(previous + extension, star_bound))
    return tuple(values)


def rank_five_facet_bound(q: int) -> int:
    """Return the first rank-five facet count allowed by the flag-UBT test."""
    if isinstance(q, bool) or not isinstance(q, int) or q < 5:
        raise ValueError("q must be an integer at least 5")
    candidate = q + 2
    while (
        q * ceil_div(candidate * (q + 1), 2)
        > 6 * (candidate * candidate - 6 * candidate + 10)
    ):
        candidate += 1
    return candidate
