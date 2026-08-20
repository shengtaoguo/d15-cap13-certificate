"""Reusable exact flag-arithmetic routines for the cap-13 certificate checker.

The module contains only the mathematical arithmetic shared by the verifier:
flag functionals, cd-coordinate conversion, trusted nonnegative factors,
Kalai convolution, generalized Dehn--Sommerville relations, and self-tests.
It has no command-line interface and contains no historical certificate logic.
"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass, field
from functools import lru_cache
from typing import Iterable, Mapping


FlagFunctional = tuple[tuple[int, int], ...]
CdCoordinates = tuple[tuple[str, int], ...]


class VerificationError(RuntimeError):
    """Raised when any fail-closed certificate gate is violated."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def normalized(values: Mapping[int, int] | Iterable[tuple[int, int]]) -> FlagFunctional:
    accumulator: dict[int, int] = {}
    items = values.items() if isinstance(values, Mapping) else values
    for raw_mask, raw_value in items:
        mask = int(raw_mask)
        value = int(raw_value)
        require(mask >= 0, "negative flag mask")
        if value:
            merged = accumulator.get(mask, 0) + value
            if merged:
                accumulator[mask] = merged
            else:
                accumulator.pop(mask, None)
    return tuple(sorted(accumulator.items()))


ONE: FlagFunctional = ((0, 1),)


def add_functionals(
    left: FlagFunctional,
    right: FlagFunctional,
    left_scale: int = 1,
    right_scale: int = 1,
) -> FlagFunctional:
    values: dict[int, int] = {}
    for mask, coefficient in left:
        value = left_scale * coefficient
        if value:
            values[mask] = value
    for mask, coefficient in right:
        value = values.get(mask, 0) + right_scale * coefficient
        if value:
            values[mask] = value
        else:
            values.pop(mask, None)
    return tuple(sorted(values.items()))


def evaluate(functional: FlagFunctional, flag_vector: Mapping[int, int]) -> int:
    return sum(coefficient * flag_vector[mask] for mask, coefficient in functional)


def dual_functional(functional: FlagFunctional, dimension: int) -> FlagFunctional:
    values: dict[int, int] = {}
    for mask, coefficient in functional:
        require(mask < (1 << dimension), "dual input mask exceeds its dimension")
        reversed_mask = 0
        for rank in range(dimension):
            if mask & (1 << rank):
                reversed_mask |= 1 << (dimension - 1 - rank)
        values[reversed_mask] = values.get(reversed_mask, 0) + coefficient
    return normalized(values)


def convolve_flag_functionals(
    left: FlagFunctional,
    left_dimension: int,
    right: FlagFunctional,
    right_dimension: int,
) -> FlagFunctional:
    """Expand Kalai convolution directly as face-chain concatenation."""
    require(left_dimension >= 0 and right_dimension >= 0, "negative factor dimension")
    left_grade = left_dimension + 1
    separator = 1 << left_dimension
    values: dict[int, int] = {}
    for left_mask, left_value in left:
        require(left_mask < (1 << left_dimension), "left mask exceeds factor dimension")
        for right_mask, right_value in right:
            require(right_mask < (1 << right_dimension), "right mask exceeds factor dimension")
            mask = left_mask | separator | (right_mask << left_grade)
            values[mask] = values.get(mask, 0) + left_value * right_value
    return normalized(values)


def add_scaled_integer(
    accumulator: dict[int, int],
    functional: FlagFunctional,
    multiplier: int,
) -> None:
    require(type(multiplier) is int, "noninteger functional multiplier")
    if multiplier == 0:
        return
    for mask, coefficient in functional:
        value = accumulator.get(mask, 0) + multiplier * coefficient
        if value:
            accumulator[mask] = value
        else:
            accumulator.pop(mask, None)


# ---------------------------------------------------------------------------
# cd words, explicit ab expansion, and independent coefficient extraction


@lru_cache(maxsize=None)
def cd_words(degree: int) -> tuple[str, ...]:
    require(degree >= 0, "negative cd degree")
    if degree == 0:
        return ("",)
    words: list[str] = []
    words.extend(word + "c" for word in cd_words(degree - 1))
    if degree >= 2:
        words.extend(word + "d" for word in cd_words(degree - 2))
    return tuple(words)


def word_degree(word: str) -> int:
    require(all(letter in "cd" for letter in word), f"invalid cd word {word!r}")
    return sum(1 if letter == "c" else 2 for letter in word)


@lru_cache(maxsize=None)
def ab_expansion(word: str) -> tuple[tuple[int, int], ...]:
    """Expand c=a+b and d=ab+ba, encoding b-positions as a mask."""
    position = 0
    terms: dict[int, int] = {0: 1}
    for letter in word:
        updated: dict[int, int] = {}
        if letter == "c":
            options = (0, 1 << position)
            width = 1
        elif letter == "d":
            options = (1 << (position + 1), 1 << position)
            width = 2
        else:
            raise VerificationError(f"invalid cd letter {letter!r}")
        for mask, coefficient in terms.items():
            for option in options:
                target = mask | option
                updated[target] = updated.get(target, 0) + coefficient
        terms = updated
        position += width
    return tuple(sorted(terms.items()))


def canonical_ab_mask(word: str) -> int:
    """Use c->a and d->ba for the declared unit-triangular pivots."""
    position = 0
    mask = 0
    for letter in word:
        if letter == "c":
            position += 1
        elif letter == "d":
            mask |= 1 << position
            position += 2
        else:
            raise VerificationError(f"invalid cd letter {letter!r}")
    return mask


@lru_cache(maxsize=None)
def flag_h_functional(mask: int) -> FlagFunctional:
    terms: list[tuple[int, int]] = []
    target_size = mask.bit_count()
    subset = mask
    while True:
        sign = -1 if ((target_size - subset.bit_count()) & 1) else 1
        terms.append((subset, sign))
        if subset == 0:
            break
        subset = (subset - 1) & mask
    return tuple(sorted(terms))


@lru_cache(maxsize=None)
def simplex_flag_vector(dimension: int) -> dict[int, int]:
    require(dimension >= 0, "negative simplex dimension")
    result = {0: 1}
    for mask in range(1, 1 << dimension):
        ranks = [rank for rank in range(dimension) if mask & (1 << rank)]
        top_size = ranks[-1] + 1
        count = math.comb(dimension + 1, top_size)
        containing_size = top_size
        for rank in reversed(ranks[:-1]):
            face_size = rank + 1
            count *= math.comb(containing_size, face_size)
            containing_size = face_size
        result[mask] = count
    return result


@dataclass
class CdBasis:
    degree: int
    words: tuple[str, ...]
    word_index: dict[str, int]
    pivots: tuple[int, ...]
    inverse_rows: tuple[tuple[tuple[int, int], ...], ...]
    simplex_coefficients: tuple[int, ...]
    extractor_cache: dict[int, FlagFunctional] = field(default_factory=dict)

    def extractor(self, index: int) -> FlagFunctional:
        require(0 <= index < len(self.words), "cd extractor index out of range")
        cached = self.extractor_cache.get(index)
        if cached is not None:
            return cached
        accumulator: FlagFunctional = ()
        for pivot_index, coefficient in self.inverse_rows[index]:
            accumulator = add_functionals(
                accumulator,
                flag_h_functional(self.pivots[pivot_index]),
                1,
                coefficient,
            )
        self.extractor_cache[index] = accumulator
        return accumulator

    def extractor_for_word(self, word: str) -> FlagFunctional:
        return self.extractor(self.word_index[word])


@lru_cache(maxsize=None)
def cd_basis(degree: int) -> CdBasis:
    words = cd_words(degree)
    pivots = tuple(canonical_ab_mask(word) for word in words)
    require(len(set(pivots)) == len(words), "duplicate canonical ab pivot")

    expansions = [dict(ab_expansion(word)) for word in words]
    matrix_rows: list[list[tuple[int, int]]] = []
    for row_index, pivot in enumerate(pivots):
        entries = [
            (column_index, expansion[pivot])
            for column_index, expansion in enumerate(expansions)
            if expansion.get(pivot, 0)
        ]
        require(
            all(column <= row_index for column, _ in entries),
            "declared cd pivot matrix is not lower triangular",
        )
        diagonal = dict(entries).get(row_index, 0)
        require(diagonal == 1, "cd pivot matrix does not have unit diagonal")
        matrix_rows.append(entries)

    # Invert the explicit unit-lower-triangular matrix in coefficient space.
    # This is intentionally separated from flag-functional construction.
    inverse: list[dict[int, int]] = []
    for row_index, entries in enumerate(matrix_rows):
        row: dict[int, int] = {row_index: 1}
        for prior_index, coefficient in entries:
            if prior_index == row_index:
                continue
            for column, prior_value in inverse[prior_index].items():
                merged = row.get(column, 0) - coefficient * prior_value
                if merged:
                    row[column] = merged
                else:
                    row.pop(column, None)
        inverse.append(row)

    # Verify M M^{-1}=I without using the flag-functional code.
    for row_index, entries in enumerate(matrix_rows):
        product: dict[int, int] = {}
        for middle, coefficient in entries:
            for column, inverse_value in inverse[middle].items():
                product[column] = product.get(column, 0) + coefficient * inverse_value
        product = {column: value for column, value in product.items() if value}
        require(product == {row_index: 1}, "exact cd pivot-matrix inversion failed")

    simplex = simplex_flag_vector(degree)
    pivot_h_values = [evaluate(flag_h_functional(mask), simplex) for mask in pivots]
    simplex_coefficients = tuple(
        sum(coefficient * pivot_h_values[column] for column, coefficient in row.items())
        for row in inverse
    )
    return CdBasis(
        degree=degree,
        words=words,
        word_index={word: index for index, word in enumerate(words)},
        pivots=pivots,
        inverse_rows=tuple(tuple(sorted(row.items())) for row in inverse),
        simplex_coefficients=simplex_coefficients,
    )


def cd_column_flag_value(word: str, flag_mask: int) -> int:
    """The f_S coordinate of the flag vector represented by one cd word."""
    return sum(
        coefficient
        for ab_mask, coefficient in ab_expansion(word)
        if not (ab_mask & ~flag_mask)
    )


def functional_cd_coordinates(
    functional: FlagFunctional, degree: int
) -> dict[str, int]:
    coordinates: dict[str, int] = {}
    for word in cd_words(degree):
        value = sum(
            coefficient * cd_column_flag_value(word, mask)
            for mask, coefficient in functional
        )
        if value:
            coordinates[word] = value
    return coordinates


# ---------------------------------------------------------------------------
# Toric g from its flag recursion and the independent cd formula


def lift_functional(functional: FlagFunctional, rank: int) -> FlagFunctional:
    bit = 1 << rank
    return tuple((mask | bit, coefficient) for mask, coefficient in functional)


@lru_cache(maxsize=None)
def toric_g_flag_functionals(dimension: int) -> tuple[FlagFunctional, ...]:
    """Compute toric g using the defining recursion over face ranks."""
    if dimension == -1:
        return (ONE,)
    require(dimension >= 0, "toric g requested below dimension -1")
    polynomial: list[FlagFunctional] = [() for _ in range(dimension + 1)]
    for face_rank in range(-1, dimension):
        lower_g = (
            (ONE,)
            if face_rank == -1
            else tuple(
                lift_functional(functional, face_rank)
                for functional in toric_g_flag_functionals(face_rank)
            )
        )
        exponent = dimension - 1 - face_rank
        coefficients = [
            math.comb(exponent, power) * ((-1) ** (exponent - power))
            for power in range(exponent + 1)
        ]
        for g_index, functional in enumerate(lower_g):
            for power, coefficient in enumerate(coefficients):
                target = g_index + power
                polynomial[target] = add_functionals(
                    polynomial[target], functional, 1, coefficient
                )
    h_rows = list(reversed(polynomial))
    require(h_rows[0] == ONE, f"toric h_0 is not one in dimension {dimension}")
    g_rows = [h_rows[0]]
    for index in range(1, dimension // 2 + 1):
        g_rows.append(add_functionals(h_rows[index], h_rows[index - 1], 1, -1))
    return tuple(g_rows)


def polynomial_product(
    left: Mapping[int, int], right: Mapping[int, int]
) -> dict[int, int]:
    result: dict[int, int] = {}
    for left_degree, left_value in left.items():
        for right_degree, right_value in right.items():
            degree = left_degree + right_degree
            result[degree] = result.get(degree, 0) + left_value * right_value
    return {degree: value for degree, value in result.items() if value}


def ballot_difference(n: int, k: int) -> int:
    first = math.comb(n, k) if 0 <= k <= n else 0
    second = math.comb(n, k - 1) if 0 <= k - 1 <= n else 0
    return first - second


def q_polynomial(n: int) -> dict[int, int]:
    return {
        degree: ((-1) ** degree) * ballot_difference(n - 1, degree)
        for degree in range((n - 1) // 2 + 1)
    }


def t_polynomial(n: int) -> dict[int, int]:
    if n % 2 == 0:
        return {}
    degree = (n - 1) // 2
    return {degree: ((-1) ** degree) * ballot_difference(n - 1, degree)}


def toric_g_contribution(word: str) -> dict[int, int]:
    c_runs: list[int] = []
    run = 0
    for letter in word:
        if letter == "c":
            run += 1
        elif letter == "d":
            c_runs.append(run)
            run = 0
        else:
            raise VerificationError(f"invalid cd word {word!r}")
    c_runs.append(run)
    d_count = len(c_runs) - 1
    polynomial = q_polynomial(c_runs[-1] + 1)
    for earlier_run in c_runs[:-1]:
        polynomial = polynomial_product(polynomial, t_polynomial(earlier_run + 1))
        if not polynomial:
            return {}
    return {degree + d_count: value for degree, value in polynomial.items()}


@lru_cache(maxsize=None)
def toric_g_cd_rows(degree: int) -> tuple[dict[str, int], ...]:
    rows: list[dict[str, int]] = [dict() for _ in range(degree // 2 + 1)]
    for word in cd_words(degree):
        for index, value in toric_g_contribution(word).items():
            rows[index][word] = value
    return tuple(rows)


# ---------------------------------------------------------------------------
# Trusted factor catalogue


def canonical_cd_coordinates(values: Mapping[str, int]) -> CdCoordinates:
    return tuple(sorted((word, int(value)) for word, value in values.items() if value))


def reverse_cd_coordinates(values: Mapping[str, int]) -> dict[str, int]:
    return {word[::-1]: value for word, value in values.items()}


@dataclass(frozen=True)
class FactorDescriptor:
    dimension: int
    kind: str
    parameter: int | str
    dual: bool
    nontrivial: int
    name: str
    cd_coordinates: CdCoordinates


@lru_cache(maxsize=None)
def factor_functional(descriptor: FactorDescriptor) -> FlagFunctional:
    dimension = descriptor.dimension
    if descriptor.kind == "toric-g":
        base = toric_g_flag_functionals(dimension)[int(descriptor.parameter)]
    elif descriptor.kind == "simplex-minimal-cd":
        basis = cd_basis(dimension)
        word = str(descriptor.parameter)
        index = basis.word_index[word]
        base = add_functionals(
            basis.extractor(index),
            ONE,
            1,
            -basis.simplex_coefficients[index],
        )
    else:
        raise VerificationError(f"unknown trusted factor kind {descriptor.kind!r}")
    return dual_functional(base, dimension) if descriptor.dual else base


@lru_cache(maxsize=None)
def factor_catalogue(dimension: int) -> tuple[FactorDescriptor, ...]:
    require(dimension >= 0, "negative factor dimension")
    basis = cd_basis(dimension)
    factors: list[FactorDescriptor] = []
    by_coordinates: dict[CdCoordinates, list[int]] = defaultdict(list)

    def push(candidate: FactorDescriptor) -> None:
        if not candidate.cd_coordinates:
            return
        collisions = by_coordinates.get(candidate.cd_coordinates, [])
        if collisions:
            functional = factor_functional(candidate)
            for prior_index in collisions:
                if factor_functional(factors[prior_index]) == functional:
                    return
        by_coordinates[candidate.cd_coordinates].append(len(factors))
        factors.append(candidate)

    for index, coordinates in enumerate(toric_g_cd_rows(dimension)):
        push(
            FactorDescriptor(
                dimension,
                "toric-g",
                index,
                False,
                index,
                f"g{index}^{dimension}",
                canonical_cd_coordinates(coordinates),
            )
        )
        push(
            FactorDescriptor(
                dimension,
                "toric-g",
                index,
                True,
                index,
                f"gd{index}^{dimension}",
                canonical_cd_coordinates(reverse_cd_coordinates(coordinates)),
            )
        )

    normalization_word = "c" * dimension
    be_index = 0
    for word, simplex_coefficient in zip(basis.words, basis.simplex_coefficients):
        coordinates = {word: 1}
        coordinates[normalization_word] = (
            coordinates.get(normalization_word, 0) - simplex_coefficient
        )
        coordinates = {key: value for key, value in coordinates.items() if value}
        if not coordinates:
            continue
        push(
            FactorDescriptor(
                dimension,
                "simplex-minimal-cd",
                word,
                False,
                1,
                f"be{be_index}^{dimension}",
                canonical_cd_coordinates(coordinates),
            )
        )
        push(
            FactorDescriptor(
                dimension,
                "simplex-minimal-cd",
                word,
                True,
                1,
                f"bed{be_index}^{dimension}",
                canonical_cd_coordinates(reverse_cd_coordinates(coordinates)),
            )
        )
        be_index += 1

    return tuple(factors)


@lru_cache(maxsize=None)
def _factor_name_index(maximum_dimension: int) -> dict[str, tuple[int, int]]:
    require(maximum_dimension >= 0, "negative maximum factor dimension")
    result: dict[str, tuple[int, int]] = {}
    for dimension in range(maximum_dimension + 1):
        for index, descriptor in enumerate(factor_catalogue(dimension)):
            require(
                descriptor.name not in result,
                f"duplicate factor name {descriptor.name!r}",
            )
            result[descriptor.name] = (dimension, index)
    return result


def factor_address(name: str, maximum_dimension: int) -> tuple[int, int]:
    """Resolve a stable factor name to its internal catalogue address."""
    require(type(name) is str and bool(name), "invalid factor name")
    index = _factor_name_index(maximum_dimension)
    require(name in index, f"unknown factor name {name!r}")
    return index[name]


def add_cd_value(target: dict[str, int], word: str, value: int) -> None:
    if not value:
        return
    merged = target.get(word, 0) + value
    if merged:
        target[word] = merged
    else:
        target.pop(word, None)


def convolve_cd_coordinates(
    left: Mapping[str, int], right: Mapping[str, int]
) -> dict[str, int]:
    """Independent quotient check for the face-chain convolution."""
    result: dict[str, int] = {}
    for left_word, left_value in left.items():
        for right_word, right_value in right.items():
            product = left_value * right_value
            add_cd_value(result, left_word + "c" + right_word, 2 * product)
            if left_word.endswith("c"):
                add_cd_value(result, left_word[:-1] + "d" + right_word, product)
            if right_word.startswith("c"):
                add_cd_value(result, left_word + "d" + right_word[1:], product)
    return result


@dataclass(frozen=True)
class Composite:
    functional: FlagFunctional
    cd_coordinates: CdCoordinates
    dimension: int
    factors: tuple[FactorDescriptor, ...]


@lru_cache(maxsize=None)
def compose_factors(composition: tuple[tuple[int, int], ...]) -> Composite:
    require(bool(composition), "empty factor composition")
    current_functional: FlagFunctional | None = None
    current_cd: dict[str, int] | None = None
    current_dimension = -1
    descriptors: list[FactorDescriptor] = []
    for dimension, factor_index in composition:
        require(dimension >= 0, "negative factor dimension in row descriptor")
        catalogue = factor_catalogue(dimension)
        require(
            0 <= factor_index < len(catalogue),
            f"factor index {factor_index} is invalid in dimension {dimension}",
        )
        descriptor = catalogue[factor_index]
        descriptors.append(descriptor)
        factor = factor_functional(descriptor)
        factor_cd = dict(descriptor.cd_coordinates)
        if current_functional is None:
            current_functional = factor
            current_cd = factor_cd
            current_dimension = dimension
        else:
            current_functional = convolve_flag_functionals(
                current_functional, current_dimension, factor, dimension
            )
            assert current_cd is not None
            current_cd = convolve_cd_coordinates(current_cd, factor_cd)
            current_dimension += dimension + 1
    assert current_functional is not None and current_cd is not None
    return Composite(
        current_functional,
        canonical_cd_coordinates(current_cd),
        current_dimension,
        tuple(descriptors),
    )


# ---------------------------------------------------------------------------
# Generalized Dehn--Sommerville decoding


def ds_relation_from_pivot(dimension: int, pivot_mask: int) -> FlagFunctional:
    require(0 <= pivot_mask < (1 << dimension), "DS pivot lies outside flag coordinates")
    choices = [
        rank
        for rank in range(dimension)
        if pivot_mask & (1 << rank)
        and (rank == dimension - 1 or pivot_mask & (1 << (rank + 1)))
    ]
    require(bool(choices), f"DS pivot {pivot_mask} is a canonical cd mask")
    removed_rank = max(choices)
    subset_mask = pivot_mask ^ (1 << removed_rank)
    marks = [-1]
    marks.extend(rank for rank in range(dimension) if subset_mask & (1 << rank))
    marks.append(dimension)
    left = max(mark for mark in marks if mark < removed_rank)
    right = min(mark for mark in marks if mark > removed_rank)
    require(
        right == removed_rank + 1 or (removed_rank == dimension - 1 and right == dimension),
        "DS pivot decoder selected a nonterminal gap rank",
    )
    relation: dict[int, int] = {}
    for rank in range(left + 1, right):
        mask = subset_mask | (1 << rank)
        relation[mask] = relation.get(mask, 0) + ((-1) ** (rank - left - 1))
    correction = 1 - ((-1) ** (right - left - 1))
    relation[subset_mask] = relation.get(subset_mask, 0) - correction
    functional = normalized(relation)
    coefficients = dict(functional)
    require(abs(coefficients.get(pivot_mask, 0)) == 1, "DS pivot coefficient is not a unit")
    require(
        all(mask < pivot_mask for mask, _ in functional if mask != pivot_mask),
        "DS relation is not lower triangular at its recorded pivot",
    )
    return functional


# ---------------------------------------------------------------------------
# Tests and full verification


def run_self_tests(
    max_degree: int = 7, exhaustive_degree: int = 7
) -> dict[str, int]:
    """Cross-check every degree, using exhaustive tests where they stay small."""
    require(max_degree >= 0, "negative self-test maximum degree")
    require(
        0 <= exhaustive_degree <= max_degree,
        "invalid exhaustive self-test degree",
    )

    def sample_indices(length: int) -> tuple[int, ...]:
        require(length > 0, "cannot sample an empty self-test collection")
        return tuple(sorted({0, length // 2, length - 1}))

    basis_pairs = 0
    toric_rows = 0
    toric_coordinate_checks = 0
    dual_rows = 0
    convolution_pairs = 0
    for degree in range(max_degree + 1):
        basis = cd_basis(degree)
        extractor_indices = (
            tuple(range(len(basis.words)))
            if degree <= exhaustive_degree
            else sample_indices(len(basis.words))
        )
        word_indices = (
            tuple(range(len(basis.words)))
            if degree <= exhaustive_degree
            else sample_indices(len(basis.words))
        )
        for extractor_index in extractor_indices:
            extractor = basis.extractor(extractor_index)
            for word_index in word_indices:
                word = basis.words[word_index]
                value = sum(
                    coefficient * cd_column_flag_value(word, mask)
                    for mask, coefficient in extractor
                )
                require(
                    value == (1 if extractor_index == word_index else 0),
                    f"cd extractor duality failed in degree {degree}",
                )
                basis_pairs += 1

        expected_g = toric_g_cd_rows(degree)
        actual_g = toric_g_flag_functionals(degree)
        require(len(expected_g) == len(actual_g), "toric g row-count mismatch")
        for expected, functional in zip(expected_g, actual_g):
            if degree <= exhaustive_degree:
                require(
                    functional_cd_coordinates(functional, degree) == expected,
                    f"toric g cd cross-check failed in degree {degree}",
                )
                toric_coordinate_checks += len(basis.words)
            else:
                for word_index in sample_indices(len(basis.words)):
                    word = basis.words[word_index]
                    actual = sum(
                        coefficient * cd_column_flag_value(word, mask)
                        for mask, coefficient in functional
                    )
                    require(
                        actual == expected.get(word, 0),
                        f"toric g cd spot-check failed in degree {degree}",
                    )
                    toric_coordinate_checks += 1
            toric_rows += 1
        catalogue = factor_catalogue(degree)
        if degree <= exhaustive_degree:
            tested_descriptors = catalogue
        else:
            selected = set(sample_indices(len(catalogue)))
            groups: dict[tuple[str, bool], list[int]] = defaultdict(list)
            for index, descriptor in enumerate(catalogue):
                groups[(descriptor.kind, descriptor.dual)].append(index)
            for indices in groups.values():
                selected.add(indices[0])
                selected.add(indices[-1])
            tested_descriptors = tuple(catalogue[index] for index in sorted(selected))
        for descriptor in tested_descriptors:
            functional = factor_functional(descriptor)
            require(
                dual_functional(dual_functional(functional, degree), degree)
                == functional,
                "duality is not an involution",
            )
            dual_rows += 1

    # Compare original-coordinate and cd-coordinate convolution on basis rows.
    for total in range(1, max_degree + 1):
        for left_degree in range(total):
            right_degree = total - left_degree - 1
            left_basis = cd_basis(left_degree)
            right_basis = cd_basis(right_degree)
            left_words = (
                left_basis.words
                if total <= exhaustive_degree
                else tuple(
                    left_basis.words[index]
                    for index in sample_indices(len(left_basis.words))
                )
            )
            right_words = (
                right_basis.words
                if total <= exhaustive_degree
                else tuple(
                    right_basis.words[index]
                    for index in sample_indices(len(right_basis.words))
                )
            )
            for left_word in left_words:
                for right_word in right_words:
                    original = convolve_flag_functionals(
                        left_basis.extractor_for_word(left_word),
                        left_degree,
                        right_basis.extractor_for_word(right_word),
                        right_degree,
                    )
                    actual = functional_cd_coordinates(original, total)
                    expected = convolve_cd_coordinates(
                        {left_word: 1}, {right_word: 1}
                    )
                    require(actual == expected, "flag/cd convolution cross-check failed")
                    convolution_pairs += 1

    return {
        "cd_extractor_pairs": basis_pairs,
        "toric_g_rows": toric_rows,
        "toric_g_coordinate_checks": toric_coordinate_checks,
        "duality_rows": dual_rows,
        "convolution_basis_pairs": convolution_pairs,
        "exhaustive_through_degree": exhaustive_degree,
        "covered_through_degree": max_degree,
    }
