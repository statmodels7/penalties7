# The Stable Log Pseudo-Determinant and Its Derivative Blocks

Evaluates \\\log\mathrm{pdet}\sum_k\lambda_kP_k\\ with
\\\operatorname{tr}(S^{+}P_k)\\ and
\\\operatorname{tr}(S^{+}P_kS^{+}P_l)\\ by the similarity transformation
[`additive_sum()`](https://statmodels7.github.io/penalties7/reference/additive_sum.md)
describes, recursing so that any number of components is served.

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

  The gap at which the components are split,
  [`additive_tol()`](https://statmodels7.github.io/penalties7/reference/additive_tol.md)
  by default.

- r_tol:

  The relative eigenvalue tolerance a group's rank is counted at,
  matching
  [`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)'s
  own `tol`.

## Value

A list of `value` (a single number), `dlog` (a numeric vector), `d2log`
(a square matrix) and `depth` (how many levels the partition needed, a
diagnostic).

## Details

A rank is read from the components normalized one by one and never from
the weighted sum, which is the convention
[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
already follows for the rank of the whole: the null space of a sum of
positive semidefinite matrices is the intersection of theirs, so a
group's rank is a property of that group and its determination is well
conditioned. The rank left to the next level is the current one less
what the dominant group peels, so it is obtained by subtraction and
never counted again.

## See also

[`additive_sum()`](https://statmodels7.github.io/penalties7/reference/additive_sum.md),
which calls it and checks its answer.
