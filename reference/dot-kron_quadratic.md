# A Quadratic Penalty Repeated Blockwise, Without Forming It

Builds the penalty of \\I_m \otimes P\\ from \\P\\ alone, and returns a
[`QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/QuadraticPenalty.md)
whose stored matrix is sparse.
[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
calls it when `blocks` is above one.

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
are `dgCMatrix` objects. `P` is checked to be positive semidefinite and
not zero by
[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
before the call.

## Details

Three quantities of a quadratic penalty come from its matrix: the rank,
the log pseudo-determinant and a basis of the null space. All three
follow from \\P\\, because the eigenvalues of \\I_m \otimes P\\ are
those of \\P\\ repeated \\m\\ times: the rank is \\mr\\, the log
pseudo-determinant is \\m\log\mathrm{pdet}(P)\\, and the null space is
\\I_m \otimes N\\. One eigendecomposition of \\P\\ supplies the
eigenvalues and the eigenvectors, and one mask selects both the rank and
the null basis. The large matrix is assembled in sparse form with
[`Matrix::kronecker()`](https://rdrr.io/pkg/Matrix/man/kronecker-methods.html)
and stored, and it is never decomposed, so the cost of the construction
is that of decomposing \\P\\.

The same identity underlies
[`parameters7::kron_identity()`](https://statmodels7.github.io/parameters7/reference/kron_identity.html),
for the covariance of grouped random effects.

## See also

[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md),
[`parameters7::kron_identity()`](https://statmodels7.github.io/parameters7/reference/kron_identity.html)
