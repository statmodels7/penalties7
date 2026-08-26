# A Quadratic Penalty Repeated Blockwise, Without Forming It

Builds the penalty of \\I_m \otimes P\\ from \\P\\ alone, and returns a
[`QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/QuadraticPenalty.md)
whose stored matrix is sparse. Called by
[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
when `blocks` is above one.

## Usage

``` r
.kron_quadratic(P, m, link_lambda, tol)
```

## Arguments

- P:

  The symmetric matrix of one block, already symmetrized.

- m:

  How many blocks, a whole number of at least two in practice.

- link_lambda:

  The linkfunctions7 link for `lambda`.

- tol:

  The relative eigenvalue tolerance for counting a zero.

## Value

A
[`QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/QuadraticPenalty.md)
with `n_coef` equal to `m * nrow(P)`, whose `P`, `DPD` and `null_basis`
are `dgCMatrix` objects. Rejects a `P` whose eigenvalues are all below
the tolerance.

## Details

A quadratic penalty needs three things from its matrix: the rank, the
log pseudo-determinant and a basis of the null space. All three follow
from \\P\\, because the eigenvalues of \\I_m \otimes P\\ are \\P\\'s
repeated \\m\\ times: the rank is \\mr\\, the log pseudo-determinant is
\\m\log\mathrm{pdet}(P)\\, and the null space is \\I_m \otimes N\\. The
large matrix is assembled with
[`Matrix::kronecker()`](https://rdrr.io/pkg/Matrix/man/kronecker-methods.html)
and stored, never decomposed.

The saving is the whole of the construction. On second differences over
ten coefficients at \\m = 200\\, the assembled route spends 4.79 s in
[`eigen()`](https://rdrr.io/r/base/eigen.html) against 2.9e-05 s here,
and 5.82 s in all against 0.001 s. The stored matrix is 0.11 MB against
30.5 MB, at a density of 0.0022, which follows from the sparse storage.

The same identity carries
[`parameters7::kron_identity()`](https://statmodels7.github.io/parameters7/reference/kron_identity.html)
on the other side of the toolkit, for the covariance of grouped random
effects.

## See also

[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md),
[`parameters7::kron_identity()`](https://statmodels7.github.io/parameters7/reference/kron_identity.html)
