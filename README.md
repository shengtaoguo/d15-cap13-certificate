# Certificate for the dimension-15 thirteen-facet theorem

This repository contains the exact certificate for the following result:

> Every convex polytope of dimension at least 15 has a three-dimensional
> face with at most 13 facets.

The mathematical argument shows that certain flag functionals are
nonnegative if every three-face has at least 14 facets. The certificate gives
an exact positive linear combination of these functionals that equals
$-1$, modulo the generalized Dehn--Sommerville relations. This is the
desired contradiction.

## Verification

Python 3.10 or newer is sufficient; no external package is needed. Verification
is single-process, is expected to use less than 2 GiB of memory, and does not
require TeX or paper-build software. Run

```sh
python3 verification/verify.py
```

from the repository root. A successful run prints a line beginning
`VERIFIED:` and exits with status `0`.

For an optional integrity check of the complete repository snapshot, run

```sh
python3 verification/verify_checksums.py
```

## Files

```text
certificate/
├── artifacts/
│   └── certificates/
│       └── d15-cap13-certificate.json   # exact certificate
├── verification/
│   ├── verify.py                        # command-line entry point
│   ├── verify_certificate.py            # certificate checker
│   ├── flag_arithmetic.py               # exact flag and cd arithmetic
│   ├── bounds.py                        # vertex and facet bounds
│   ├── verify_checksums.py              # optional integrity checker
│   └── verification-report.json         # accepted verification result
├── repository-manifest.json             # claim and artifact metadata
└── SHA256SUMS                           # repository file hashes
```

Only `verification/verify.py` is intended to be run directly. The verification
flow is

```text
verify.py -> verify_certificate.py -> {flag_arithmetic.py, bounds.py}
                    |
                    v
       d15-cap13-certificate.json
```

## Certificate identity

The JSON file records the identity

$$
1 + \sum_i \lambda_i L_i + \sum_j \mu_j D_j = 0.
$$

Here the $L_i$ are valid nonnegative flag-functionals used in the proof,
every $\lambda_i$ is a strictly positive rational number, and the $D_j$ are
generalized Dehn--Sommerville relations. The certificate contains 987 rows
$L_i$ and 13,667 nonzero coefficients $\mu_j$. The checker reconstructs
the identity in all $2^{15}=32,768$ flag coordinates and verifies that
every coordinate is exactly zero.

To avoid repeating large rational denominators, the JSON file stores the
equivalent primitive integer identity

$$
N + \sum_i a_i L_i + \sum_j b_j D_j = 0,
$$

where $N>0$, $\lambda_i=a_i/N$, and $\mu_j=b_j/N$. Thus every $a_i$
is strictly positive, while the $b_j$ may have either sign. The checker
uses only integer arithmetic.

## Certificate format

The file has the following hierarchy:

```text
certificate
├── schema
├── parameters
├── normalization                         N
├── inequality_terms                      the 987 terms a_i L_i
│   ├── row
│   │   ├── anchor                        omitted when unconditional
│   │   │   ├── family                    type of anchor inequality
│   │   │   └── dimension                 dimension of the anchor
│   │   └── factors                       ordered factor names
│   └── coefficient                       positive integer a_i
└── dehn_sommerville_terms                the 13,667 terms b_j D_j
    ├── pivot_mask                        encoding of its leading flag coordinate
    └── coefficient                       nonzero integer b_j
```

The top-level fields are:

| Field | Meaning |
| --- | --- |
| `schema` | The version of the certificate format. |
| `parameters.dimension` | The ambient polytope dimension, here 15. |
| `parameters.face_dimension` | The dimension of the faces in the theorem, here 3. |
| `parameters.minimum_facets` | The contradiction assumption: every three-face has at least 14 facets. |
| `normalization` | The positive common denominator $N$. It is stored as a 517-digit decimal string. |
| `inequality_terms` | The 987 positive terms $a_iL_i$. |
| `dehn_sommerville_terms` | The 13,667 nonzero terms $b_jD_j$. |

The large integers are written as JSON strings rather than JSON numbers. This
preserves them exactly in software whose numeric type cannot represent
500-digit integers.

### Inequality terms

An entry of `inequality_terms` represents one summand $a_iL_i$:

```json
{
  "row": {
    "anchor": {
      "family": "CAP",
      "dimension": 3
    },
    "factors": ["g1^2", "g1^2", "g0^0", "g0^0", "g0^0", "g0^0", "g0^0", "g0^0"]
  },
  "coefficient": "..."
}
```

The fields mean:

| Field | Meaning |
| --- | --- |
| `row.anchor` | The conditional anchor. It is omitted for an unconditional row. |
| `row.anchor.family` | The family of the anchor inequality. |
| `row.anchor.dimension` | The dimension of the anchor factor. |
| `row.factors` | Stable names of the remaining factors, in convolution order. |
| `coefficient` | The positive integer $a_i=N\lambda_i$. |

Writing $r$ for `row.anchor.dimension`, the row types are:

| Row type or anchor family | Number | Allowed dimension | Starting inequality |
| --- | ---: | ---: | --- |
| unconditional | 748 | no anchor | The named `factors` describe the entire row. |
| `CAP` | 123 | 3 | $f_2(Q)-14\geq 0$: a three-face $Q$ has at least 14 facets. |
| `K4` | 38 | 4 | $f_3(Q)-15\geq 0$: a four-face $Q$ has at least 15 facets. |
| `KUBT` | 61 | 5--15 | $f_{r-1}(Q)-(r+19)\geq 0$: the higher-dimensional facet bounds. |
| `V` | 17 | 3--15 | $f_0(Q)-V_r\geq 0$: the vertex bound $V_r$ computed by `bounds.py`. |

Here $f_t(Q)$ is the number of $t$-dimensional faces of $Q$. The
checker derives every constant in this table; the JSON does not store those
constants separately. For the `V` family, the values for dimensions 3
through 15 are

$$
(V_3,\ldots,V_{15})
=(9,15,21,27,33,39,45,51,56,61,66,71,76).
$$

### Factors

`row.factors` is the ordered list of unconditional factors. Each factor is
recorded by its stable catalogue name rather than by a numerical list index.
The naming conventions are:

| Catalogue kind | Meaning |
| --- | --- |
| `g` | A coefficient functional of the toric $g$-polynomial. |
| `gd` | The corresponding toric $g$-functional with face dimensions reversed. |
| `be` | A simplex-minimal $cd$-coordinate functional. |
| `bed` | The corresponding simplex-minimal functional with face dimensions reversed. |

In a name such as `g1^4`, the superscript 4 is the factor dimension. Thus
`g1^4` is the dimension-four toric $g_1$-functional, while `gd1^4` is its
dual. Every name is unique across the complete dimension-0-through-15
catalogue. The checker resolves the name and regenerates the factor's exact
integer flag coordinates.

The paper supplies the mathematical justification that the catalogue factors
and the declared anchor rows are nonnegative in the required setting. The
checker does not reprove those theorems; it checks that every recorded name
and anchor belongs to the closed allowed list and then verifies the exact
certificate identity.

Order matters: the checker applies the factors from left to right in exactly
the recorded order.

If an unconditional row has $m$ factors of dimensions
$d_1,\ldots,d_m$, its final dimension is

$$
d_1+\cdots+d_m+(m-1).
$$

If a conditional row has an anchor of dimension $r$ followed by those
$m$ factors, its final dimension is

$$
r+d_1+\cdots+d_m+m.
$$

In either case, the final dimension must equal 15.

For example, the complete dimension-two catalogue is:

| Factor name | Flag-functional |
| --- | --- |
| `g0^2` | $1$ |
| `g1^2` | $f_0-3$ |
| `gd1^2` | $f_1-3$ |

Thus the name `g1^2` directly selects the nonnegative polygon functional
$g_1^2=f_0-3$.

For a complete example, consider the following row descriptor from the
certificate:

```json
{
  "anchor": {
    "family": "V",
    "dimension": 3
  },
  "factors": ["g1^2", "g1^2", "g1^4", "g0^0"]
}
```

The entries resolve as follows:

| JSON data | Selected factor |
| --- | --- |
| `"family": "V", "dimension": 3` | The anchor $f_0-9$. |
| `"g1^2"` | $g_1^2=f_0-3$. |
| `"g1^2"` | A second copy of $g_1^2=f_0-3$. |
| `"g1^4"` | $g_1^4=f_0-5$. |
| `"g0^0"` | $g_0^0=1$. |

Writing $\star$ for convolution, the checker constructs

$$
L_i=(f_0-9)\star(f_0-3)\star(f_0-3)\star(f_0-5)\star 1.
$$

Each occurrence of $f_0$ belongs to its own factor dimension. This row has
dimension

$$
3+2+2+4+0+4=15,
$$

where the final 4 counts the four convolutions.

The optional anchor and the ordered `factors` uniquely determine $L_i$.

### Dehn--Sommerville terms

For a subset $S\subseteq\{0,\ldots,14\}$, let $f_S$ be the flag number
that counts chains of proper faces whose dimensions are exactly the elements
of $S$. The 32,768 flag coordinates are indexed by these subsets.

An entry of `dehn_sommerville_terms` represents one summand $b_jD_j$:

```json
{
  "pivot_mask": 32641,
  "coefficient": "..."
}
```

The fields mean:

| Field | Meaning |
| --- | --- |
| `pivot_mask` | A bit mask encoding the leading flag coordinate of $D_j$. |
| `coefficient` | The nonzero integer $b_j=N\mu_j$, which may have either sign. |

The bit in position $t$ is 1 exactly when $t\in S$; equivalently,

$$
\mathrm{pivot\_mask}=\sum_{t\in S}2^t.
$$

Thus `32641` encodes
$S=\{0,7,8,9,10,11,12,13,14\}$. For this pivot the checker reconstructs

$$
D_j=f_{\{0,7,8,9,10,11,12,13,14\}}
    -2f_{\{0,7,8,9,10,11,12,13\}}.
$$

The `pivot_mask` is therefore not an arbitrary identifier and does not store
the relation's coefficients. Here a pivot means the leading coordinate of
the relation: its coefficient is a unit, and every other mask in the relation
is smaller. Together with dimension 15, the mask determines a unique
lower-triangular generalized Dehn--Sommerville relation, which the checker
reconstructs exactly.

## Verification checks

The checker verifies the computational layer of the argument. It uses the
paper's nonnegativity theorems as mathematical inputs rather than attempting
to prove them. It reconstructs every recorded row and then verifies:

1. the schema, parameters, primitive normalization, and certificate file hash;
2. the formulas for all conditional facet and vertex bounds;
3. the convolution order and dimension of every row;
4. the strict positivity of all 987 integer coefficients $a_i$;
5. zero residual in the 987 degree-15 `cd`-coordinates;
6. zero residual in all 32,768 original flag coordinates.

The checker also tests eight deliberately corrupted copies of the
certificate. The line `rejected: 8` means that all eight corruptions were
detected; it is part of a successful verification.

Before checking the certificate, the checker cross-checks its independent
flag- and `cd`-coordinate implementations exhaustively through degree 7. The
subsequent certificate replay independently reconstructs the degree-15
identity in both the `cd` quotient and all $2^{15}=32,768$ original flag
coordinates.

The repository verifies the finished exact certificate. It does not rerun
the numerical search that found the 987 rows, and it does not require an
optimizer or the complete row catalogue.
