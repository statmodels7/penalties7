# Methods of penalty_dhessian_beta()

The method of the base class signals an error that names the penalty; it
is the method that
[`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
and
[`mcp_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
reach. The quadratic, additive and structured branches return a zero
matrix, because their Hessian is free of the coefficients. The separable
branch returns zero when
[`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md)
is `TRUE` for it, and \\-D'\mathrm{diag}(\ell^{(yyy)}(D\beta)\odot
Dv)D\\ otherwise, with a block-diagonal middle matrix for a multivariate
parent. It signals an error for a kinked parent and for a multivariate
parent that cannot supply its third response derivative.

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
