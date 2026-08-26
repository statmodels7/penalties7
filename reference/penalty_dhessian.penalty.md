# What the Base Class Answers to the Marginal Generics

One page for the base class's answers to the four generics a marginal
criterion asks beyond the second order.
[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md),
[`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md)
and
[`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md)
each reject, naming the penalty and the generic;
[`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md)
answers `FALSE`.

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object of a class registering none of these.

- beta:

  A numeric vector of coefficients. Unused.

- theta:

  A named list of hyperparameter values. Unused.

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

A marginal criterion is a Laplace approximation and asks for derivatives
beyond the second. A penalty that does not supply them cannot be
estimated by one, and reporting that is better than a criterion
assembled from a quantity nobody wrote.
[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)'s
message says so in full; the other two name the generic alone, a caller
who has reached the second having already read the first.

[`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md)'s
`FALSE` is the conservative default: a penalty that says nothing about
its third \\\beta\\-derivative is treated as having one.

## See also

[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
and its siblings for the generics,
[`penalty_dhessian.QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.QuadraticPenalty.md)
for a branch that answers.
