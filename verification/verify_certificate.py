#!/usr/bin/env python3
"""Independent exact checker for the dimension-15 cap-13 certificate.

The checker uses only the Python standard library and the two transparent
local modules bounds and flag_arithmetic. It reads the primitive integer
certificate, rederives every q=14 anchor and suffix convolution, and verifies
the identity in both the 987-dimensional cd quotient and all 2^15 original
flag coordinates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from collections import Counter
from copy import deepcopy
from math import gcd
from pathlib import Path

import bounds
import flag_arithmetic as fa


SCHEMA = "d15-cap13-independent-validity-check-v2"
CERTIFICATE_SCHEMA = "d15-cap13-exact-certificate-v3"
EXPECTED_CERTIFICATE_SHA256 = (
    "695ef9af1b23a325b01186b7d9492448b41d9ef2a06a4fe9af5f7767000b5ff4"
)
EXPECTED_FAMILY_CENSUS = {
    "CAP:r3": 123,
    "K4:r4": 38,
    "KUBT:r5": 13,
    "KUBT:r6": 12,
    "KUBT:r7": 11,
    "KUBT:r8": 8,
    "KUBT:r9": 6,
    "KUBT:r10": 4,
    "KUBT:r11": 3,
    "KUBT:r12": 1,
    "KUBT:r13": 1,
    "KUBT:r14": 1,
    "KUBT:r15": 1,
    "V:r3": 11,
    "V:r4": 2,
    "V:r6": 1,
    "V:r7": 2,
    "V:r12": 1,
    "unconditional": 748,
}


class VerificationError(RuntimeError):
    """Raised when any fail-closed certificate gate is violated."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_integer(value: object, label: str) -> int:
    require(type(value) is str, f"{label} is not an integer string")
    try:
        result = int(value)
    except ValueError as error:
        raise VerificationError(f"malformed integer in {label}") from error
    require(str(result) == value, f"noncanonical integer string in {label}")
    return result


def parse_factor_names(raw: object, *, allow_empty: bool) -> tuple[str, ...]:
    """Parse an ordered list of stable factor names."""
    require(isinstance(raw, list), "row factors is not a list")
    if not raw:
        require(allow_empty, "unconditional factor list is empty")
        return ()
    result: list[str] = []
    for item in raw:
        require(type(item) is str and bool(item), "malformed factor name")
        result.append(item)
    return tuple(result)


def compose_named_factors(names: tuple[str, ...], maximum_dimension: int):
    addresses = tuple(
        fa.factor_address(name, maximum_dimension) for name in names
    )
    composite = fa.compose_factors(addresses)
    require(
        tuple(descriptor.name for descriptor in composite.factors) == names,
        "factor name resolution changed",
    )
    return composite


def anchor_functional(kind: str, rank: int, q: int):
    require(q == 14, "certificate anchor threshold is not q=14")
    if kind == "CAP":
        require(rank == 3, "CAP anchor has wrong rank")
        return fa.normalized({1 << 2: 1, 0: -q})
    if kind == "V":
        require(3 <= rank <= 15, "V anchor rank is out of range")
        constant = bounds.vertex_bounds(q, rank)[rank - 3]
        return fa.normalized({1 << 0: 1, 0: -constant})
    if kind == "K4":
        require(rank == 4, "K4 anchor has wrong rank")
        return fa.normalized({1 << 3: 1, 0: -(q + 1)})
    if kind == "KUBT":
        require(5 <= rank <= 15, "KUBT anchor rank is out of range")
        constant = bounds.rank_five_facet_bound(q) + rank - 5
        return fa.normalized({1 << (rank - 1): 1, 0: -constant})
    raise VerificationError(f"unrecognized T56 anchor kind {kind!r}")


def expected_anchor_cd(kind: str, rank: int, q: int) -> dict[str, int]:
    normalization = "c" * rank
    if kind == "CAP":
        return {"cd": 1, normalization: -(q - 2)}
    if kind == "V":
        constant = bounds.vertex_bounds(q, rank)[rank - 3]
        return {"d" + "c" * (rank - 2): 1, normalization: -(constant - 2)}
    if kind in {"K4", "KUBT"}:
        constant = q + 1 if kind == "K4" else bounds.rank_five_facet_bound(q) + rank - 5
        return {"c" * (rank - 2) + "d": 1, normalization: -(constant - 2)}
    raise VerificationError(f"unrecognized T56 anchor kind {kind!r}")


def decode_row(row: object, dimension: int, q: int):
    require(isinstance(row, dict), "row descriptor is not an object")
    if "anchor" not in row:
        require(set(row) == {"factors"}, "unconditional row has unexpected fields")
        factor_names = parse_factor_names(row.get("factors"), allow_empty=False)
        composite = compose_named_factors(factor_names, dimension)
        require(
            composite.dimension == dimension,
            "unconditional row has wrong dimension",
        )
        return (
            composite.functional,
            composite.cd_coordinates,
            composite.factors,
            "unconditional",
        )

    require(
        set(row) == {"anchor", "factors"},
        "conditional row has unexpected fields",
    )
    anchor_descriptor = row.get("anchor")
    require(isinstance(anchor_descriptor, dict), "anchor is not an object")
    require(
        set(anchor_descriptor) == {"family", "dimension"},
        "anchor has unexpected fields",
    )
    family = anchor_descriptor.get("family")
    require(
        family in {"CAP", "V", "K4", "KUBT"},
        f"unrecognized row family {family!r}",
    )
    anchor_dimension = anchor_descriptor.get("dimension")
    require(
        type(anchor_dimension) is int,
        "anchor dimension is not an integer",
    )
    anchor = anchor_functional(family, anchor_dimension, q)
    anchor_cd = fa.functional_cd_coordinates(anchor, anchor_dimension)
    require(
        anchor_cd == expected_anchor_cd(family, anchor_dimension, q),
        "anchor closed cd formula failed",
    )
    factor_names = parse_factor_names(row.get("factors"), allow_empty=True)
    if not factor_names:
        require(
            anchor_dimension == dimension,
            "empty suffix occurs below ambient dimension",
        )
        functional = anchor
        coordinates = anchor_cd
        factors = ()
    else:
        suffix = compose_named_factors(factor_names, dimension)
        require(
            suffix.dimension == dimension - anchor_dimension - 1,
            "conditional suffix has wrong dimension",
        )
        functional = fa.convolve_flag_functionals(
            anchor, anchor_dimension, suffix.functional, suffix.dimension
        )
        coordinates = fa.convolve_cd_coordinates(
            anchor_cd, dict(suffix.cd_coordinates)
        )
        factors = suffix.factors
    return (
        functional,
        fa.canonical_cd_coordinates(coordinates),
        factors,
        f"{family}:r{anchor_dimension}",
    )


def verify_payload(
    payload: object, *, run_self_tests: bool, verbose: bool
) -> dict[str, object]:
    started = time.perf_counter()
    require(isinstance(payload, dict), "certificate root is not an object")
    require(
        set(payload)
        == {
            "schema",
            "parameters",
            "normalization",
            "inequality_terms",
            "dehn_sommerville_terms",
        },
        "certificate root has unexpected fields",
    )
    require(payload.get("schema") == CERTIFICATE_SCHEMA, "certificate schema changed")
    require(
        payload.get("parameters")
        == {
            "dimension": 15,
            "face_dimension": 3,
            "minimum_facets": 14,
        },
        "certificate parameters changed",
    )
    normalization = parse_integer(payload.get("normalization"), "normalization")
    require(normalization > 0, "normalization is not strictly positive")

    tests = fa.run_self_tests(max_degree=7) if run_self_tests else {}
    if run_self_tests:
        for kind, ranks in {
            "CAP": (3,),
            "V": tuple(range(3, 16)),
            "K4": (4,),
            "KUBT": tuple(range(5, 16)),
        }.items():
            for rank in ranks:
                actual = fa.functional_cd_coordinates(
                    anchor_functional(kind, rank, 14), rank
                )
                require(
                    actual == expected_anchor_cd(kind, rank, 14),
                    "T56 anchor self-test failed",
                )
        tests = {**tests, "t56_anchor_cd_rows": 26}

    inequality_terms = payload.get("inequality_terms")
    ds_terms = payload.get("dehn_sommerville_terms")
    require(
        isinstance(inequality_terms, list) and len(inequality_terms) == 987,
        "inequality term count mismatch",
    )
    require(
        isinstance(ds_terms, list) and len(ds_terms) == 13667,
        "Dehn--Sommerville term count mismatch",
    )

    quotient_residual: dict[str, int] = {"c" * 15: normalization}
    original_residual: dict[int, int] = {}
    fa.add_scaled_integer(original_residual, fa.ONE, normalization)
    coefficient_gcd = normalization
    family_census: Counter[str] = Counter()
    factor_occurrences: Counter[str] = Counter()
    unique_factors: set[str] = set()

    for position, entry in enumerate(inequality_terms):
        require(
            isinstance(entry, dict)
            and set(entry) == {"coefficient", "row"},
            f"inequality_terms[{position}] is malformed",
        )
        coefficient = parse_integer(
            entry.get("coefficient"),
            f"inequality_terms[{position}].coefficient",
        )
        require(coefficient > 0, "inequality coefficient is not strictly positive")
        coefficient_gcd = gcd(coefficient_gcd, coefficient)
        functional, coordinates, factors, census_key = decode_row(
            entry.get("row"), 15, 14
        )
        family_census[census_key] += 1
        factor_names = parse_factor_names(
            entry["row"].get("factors"), allow_empty=True
        )
        for name, descriptor in zip(factor_names, factors, strict=True):
            unique_factors.add(name)
            factor_occurrences[
                descriptor.kind + ("-dual" if descriptor.dual else "")
            ] += 1
        for word, row_coefficient in coordinates:
            value = quotient_residual.get(word, 0) + coefficient * row_coefficient
            if value:
                quotient_residual[word] = value
            else:
                quotient_residual.pop(word, None)
        fa.add_scaled_integer(original_residual, functional, coefficient)
        if verbose and (position + 1) % 100 == 0:
            print(
                f"validated {position + 1}/987 T56 inequality rows",
                flush=True,
            )

    require(
        dict(sorted(family_census.items())) == EXPECTED_FAMILY_CENSUS,
        "derived family census changed",
    )
    require(not quotient_residual, "exact quotient identity has nonzero residual")

    canonical_masks = {fa.canonical_ab_mask(word) for word in fa.cd_words(15)}
    pivots: set[int] = set()
    for position, entry in enumerate(ds_terms):
        require(
            isinstance(entry, dict)
            and set(entry) == {"coefficient", "pivot_mask"},
            f"dehn_sommerville_terms[{position}] is malformed",
        )
        pivot = entry.get("pivot_mask")
        require(type(pivot) is int and 0 <= pivot < (1 << 15), "invalid DS pivot")
        require(pivot not in canonical_masks, "DS support contains a canonical pivot")
        require(pivot not in pivots, "duplicate DS pivot")
        pivots.add(pivot)
        coefficient = parse_integer(
            entry.get("coefficient"),
            f"dehn_sommerville_terms[{position}].coefficient",
        )
        require(coefficient != 0, "stored DS coefficient is zero")
        coefficient_gcd = gcd(coefficient_gcd, abs(coefficient))
        fa.add_scaled_integer(
            original_residual,
            fa.ds_relation_from_pivot(15, pivot),
            coefficient,
        )
        if verbose and (position + 1) % 2000 == 0:
            print(f"validated {position + 1}/13667 DS rows", flush=True)
    require(coefficient_gcd == 1, "integer certificate is not primitively normalized")
    require(
        not original_residual,
        "full original-coordinate identity has nonzero residual",
    )

    return {
        "schema": SCHEMA,
        "status": "VERIFIED",
        "dimension": 15,
        "q": 14,
        "cap": 13,
        "original_flag_coordinates": 1 << 15,
        "quotient_coordinates": 987,
        "integer_normalization_digits": len(str(normalization)),
        "inequality_rows": len(inequality_terms),
        "strictly_positive_inequality_multipliers": len(inequality_terms),
        "dehn_sommerville_rows": len(ds_terms),
        "nonzero_dehn_sommerville_multipliers": len(ds_terms),
        "family_census": dict(sorted(family_census.items())),
        "unique_suffix_factors": len(unique_factors),
        "factor_occurrences": dict(sorted(factor_occurrences.items())),
        "self_tests": tests,
        "checks": {
            "closed_factor_catalogue": True,
            "all_anchor_formulas_rederived": True,
            "all_convolution_dimensions_valid": True,
            "all_multipliers_strictly_positive": True,
            "primitive_integer_normalization": True,
            "quotient_residual_zero": True,
            "all_ds_rows_match_declared_formula": True,
            "original_coordinate_residual_zero": True,
            "floating_point_used": False,
            "optimizer_modules_imported": False,
        },
        "elapsed_seconds": time.perf_counter() - started,
    }


def mutation_tests(payload: dict[str, object]) -> dict[str, object]:
    def first_conditional(data: dict[str, object]) -> dict[str, object]:
        return next(
            entry
            for entry in data["inequality_terms"]
            if "anchor" in entry["row"]
        )

    mutations = []

    def add(name, mutate):
        candidate = deepcopy(payload)
        mutate(candidate)
        try:
            verify_payload(candidate, run_self_tests=False, verbose=False)
        except (VerificationError, fa.VerificationError, ValueError):
            mutations.append({"name": name, "status": "REJECTED"})
            return
        raise VerificationError(f"mutation was not rejected: {name}")

    add(
        "wrong dimension",
        lambda data: data["parameters"].__setitem__("dimension", 14),
    )
    add(
        "wrong threshold",
        lambda data: data["parameters"].__setitem__("minimum_facets", 15),
    )
    add("zero normalization", lambda data: data.__setitem__("normalization", "0"))
    add(
        "negative inequality coefficient",
        lambda data: data["inequality_terms"][0].__setitem__(
            "coefficient", "-1"
        ),
    )
    add(
        "duplicate DS pivot",
        lambda data: data["dehn_sommerville_terms"][1].__setitem__(
            "pivot_mask", data["dehn_sommerville_terms"][0]["pivot_mask"]
        ),
    )
    add(
        "unknown factor name",
        lambda data: data["inequality_terms"][0]["row"]["factors"].__setitem__(
            0, "UNKNOWN"
        ),
    )
    add(
        "unknown row family",
        lambda data: first_conditional(data)["row"]["anchor"].__setitem__(
            "family", "UNKNOWN"
        ),
    )
    add(
        "zero DS coefficient",
        lambda data: data["dehn_sommerville_terms"][0].__setitem__(
            "coefficient", "0"
        ),
    )
    return {"status": "PASS", "rejected": len(mutations), "mutations": mutations}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--expected-sha256", default=EXPECTED_CERTIFICATE_SHA256)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--mutation-tests", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    certificate = args.certificate.resolve()
    certificate_bytes = certificate.read_bytes()
    observed_sha256 = hashlib.sha256(certificate_bytes).hexdigest()
    require(observed_sha256 == args.expected_sha256, "certificate SHA-256 mismatch")
    payload = json.loads(certificate_bytes.decode("utf-8"))
    report = verify_payload(payload, run_self_tests=True, verbose=args.verbose)
    report.update(
        {
            "certificate_sha256": observed_sha256,
            "checker_sha256": sha256_file(Path(__file__).resolve()),
            "arithmetic_sha256": sha256_file(Path(fa.__file__).resolve()),
            "bounds_sha256": sha256_file(Path(bounds.__file__).resolve()),
        }
    )
    if args.mutation_tests:
        report["mutation_tests"] = mutation_tests(payload)
    if args.report is not None:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print(
        "VERIFIED: 987 valid positive T56 inequalities; "
        "13,667 valid DS terms; zero exact residual in 32,768 flag coordinates"
    )
    print(json.dumps(report, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
