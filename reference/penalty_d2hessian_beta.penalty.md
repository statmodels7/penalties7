# Methods of penalty_d2hessian_beta() and penalty_dhessian_beta_theta()

The methods of the base class for both generics signal an error that
names the penalty; they are the methods that
[`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
and
[`mcp_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
reach. The quadratic, additive and structured branches
([`QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/QuadraticPenalty.md),
[`AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/AdditivePenalty.md)
and
[`StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/StructuredPenalty.md))
share one method per generic, which returns a zero matrix of side
`pen@n_coef` from
[`penalty_d2hessian_beta()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.md)
and a list of such zero matrices, one per hyperparameter, from
[`penalty_dhessian_beta_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.md),
because their Hessian is free of the coefficients. The separable branch
returns the same zeros when
[`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md)
is `TRUE` for it, and otherwise the closed forms given on the page of
[`penalty_d2hessian_beta()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.md),
contracted per block for a multivariate parent. It signals an error for
a kinked parent and for a multivariate parent that cannot supply its
fourth or mixed third response derivative.

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
for the generics,
[`penalty_dhessian_beta.penalty()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian_beta.penalty.md)
for the methods of the derivative of first order.
