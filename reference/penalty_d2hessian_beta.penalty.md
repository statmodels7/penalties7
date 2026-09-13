# What Each Branch Answers to the Second Movement of the Hessian

The base class rejects both, naming the penalty. The quadratic, additive
and structured branches return zero matrices, their Hessian being free
of the coefficients. The separable branch returns zero where its parent
is quadratic in the argument, the closed forms of
[`penalty_d2hessian_beta()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.md)
otherwise, and rejects for a kinked parent and for a multivariate parent
that is not quadratic.

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- beta:

  A numeric vector of length `pen@n_coef`.

- theta:

  A named list of hyperparameter values.

- v, w:

  Numeric vectors of length `pen@n_coef`, the directions.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A square base matrix of side `pen@n_coef`, a list of them keyed by
`pen@params`, or an error.

## See also

[`penalty_d2hessian_beta()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.md)
for the generics.
