# A Quadratic Penalty Repeated Blockwise, Without Forming It

The penalty of \\I_m \otimes P\\, built from \\P\\ alone.

## Usage

``` r
.kron_quadratic(P, m, link_lambda, tol)
```

## Arguments

- P:

  The symmetric matrix of one block.

- m:

  How many blocks.

- link_lambda:

  The link for the hyperparameter.

- tol:

  The relative tolerance for a zero eigenvalue.

## Value

A `QuadraticPenalty`.

## Details

What the constructor needs from the matrix is its rank, its log pseudo-
determinant and a basis of its null space, and all three of them follow
from \\P\\: the eigenvalues of \\I_m \otimes P\\ are \\P\\'s repeated
\\m\\ times, so the rank is \\m\\r\\, the log pseudo- determinant is
\\m\log\mathrm{pdet}(P)\\ and the null space is \\I_m \otimes N\\. The
big matrix is therefore never decomposed.

It is what one smooth per level of a factor needs, and the saving is the
whole of the construction: at \\m = 200\\ over a basis of ten, the
eigendecomposition of the assembled \\1800 \times 1800\\ matrix costs
4.50 seconds and the one of \\P\\ costs nothing measurable. The stored
matrix is sparse besides – 25.9 MB dense at a density of 0.0005 – which
follows rather than being the point.

The same identity is what
[`kron_identity`](https://statmodels7.github.io/parameters7/reference/kron_identity.html)
uses on the other side of the toolkit, for the covariance of grouped
random effects.

## See also

[`quadratic_penalty`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
