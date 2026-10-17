# Base Class Methods of the Marginal Generics

The methods of the base class cover the four generics that a marginal
criterion requires beyond the second order.
[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md),
[`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md)
and
[`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md)
each signal an error that names the penalty and the generic;
[`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md)
returns `FALSE`.
[`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
and
[`mcp_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
register none of the four and reach these methods.

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object of a class that registers none of these methods.

- beta:

  A numeric vector of coefficients. The methods do not use it.

- theta:

  A named list of hyperparameter values. The methods do not use it.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md),
[`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md)
and
[`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md)
signal an error.
[`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md)
returns the single logical `FALSE`.

## Details

A marginal criterion based on a Laplace approximation requires
derivatives of the penalty beyond the second order. For a penalty that
does not supply them, these methods signal an error, so the
hyperparameters of such a penalty cannot be estimated by a marginal
criterion. The message of
[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
states that consequence; the messages of
[`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md)
and
[`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md)
name only the penalty and the generic.

The `FALSE` of
[`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md)
is the conservative default: a penalty class without a method of its own
is treated as having a third \\\beta\\-derivative.

## See also

[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
and its siblings for the generics,
[`penalty_dhessian.QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.QuadraticPenalty.md)
for a branch that supplies them.
