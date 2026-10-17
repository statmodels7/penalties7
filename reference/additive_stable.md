# The Stable Log Pseudo-Determinant and Its Derivative Blocks

Evaluates \\\log\mathrm{pdet}\sum_k\lambda_kP_k\\ together with
\\\operatorname{tr}(S^{+}P_k)\\ and
\\\operatorname{tr}(S^{+}P_kS^{+}P_l)\\ by the similarity transformation
that
[`additive_sum()`](https://statmodels7.github.io/penalties7/reference/additive_sum.md)
describes. The function recurses on the Schur complement, so it handles
any number of components.

## Usage

``` r
additive_stable(mats, lam, p_rank, d_tol = additive_tol(), r_tol = 1e-10)
```

## Arguments

- mats:

  The component matrices, as
  [`AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/AdditivePenalty.md)
  stores them.

- lam:

  The smoothing parameters, in the same order.

- p_rank:

  The rank of the sum, fixed at construction.

- d_tol:

  The threshold at which the components are split,
  [`additive_tol()`](https://statmodels7.github.io/penalties7/reference/additive_tol.md)
  by default.

- r_tol:

  The relative eigenvalue tolerance at which the rank of a group is
  counted, `1e-10` by default.
  [`additive_sum()`](https://statmodels7.github.io/penalties7/reference/additive_sum.md)
  passes the `tol` stored in the
  [`AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/AdditivePenalty.md)
  object, which is the tolerance at which
  [`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
  counted the rank of the sum.

## Value

A list of `value` (a single number), `dlog` (a numeric vector), `d2log`
(a square matrix) and `depth` (the number of levels of the partition, a
diagnostic).

## Details

The rank of a group of components is counted on the sum of its
components after each has been divided by its largest absolute entry,
and never on the weighted sum.
[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
counts the rank of the whole penalty in the same way and at the same
tolerance. The null space of a sum of positive semidefinite matrices is
the intersection of their null spaces, so the rank of a group is a
property of that group and is counted without the spread of the
parameters. The rank passed to the next level is the current rank less
the rank of the dominant group, so it is obtained by subtraction and
never counted again.

A level returns `NaN` for the value and for both derivative blocks when
its weighted sum or reduction is not finite, when the Cholesky
factorization of \\A\\ fails, or, at the last level, when the smallest
of the eigenvalues that are kept is not positive; a `NaN` at a deeper
level is passed up.

## See also

[`additive_sum()`](https://statmodels7.github.io/penalties7/reference/additive_sum.md),
which calls it and checks its result.
