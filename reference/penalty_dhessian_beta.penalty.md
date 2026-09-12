# What Each Branch Answers to penalty_dhessian_beta()

The base class rejects, naming the penalty. The quadratic, additive and
structured branches return a zero matrix, their Hessian being free of
the coefficients. The separable branch returns zero where its parent is
quadratic in the argument and
\\-D'\mathrm{diag}(\ell^{(yyy)}(D\beta)\odot Dv)D\\ otherwise, and
rejects for a kinked parent and for a multivariate parent that is not
quadratic.

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- beta:

  A numeric vector of length `pen@n_coef`.

- theta:

  A named list of hyperparameter values.

- v:

  A numeric vector of length `pen@n_coef`, the direction.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A square base matrix of side `pen@n_coef`, or an error.

## See also

[`penalty_dhessian_beta()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian_beta.md)
for the generic.
