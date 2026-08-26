# Marginal Derivatives of a Quadratic Penalty

One page for the branch's answers to the four generics a marginal
criterion asks beyond the second order. \\S = \lambda D'PD\\ is linear
in \\\lambda\\ and free of the coefficients, so
[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
returns \\D'PD\\,
[`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md)
and
[`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md)
return zeros, and
[`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md)
returns `TRUE`.

## Arguments

- pen:

  A
  [`QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/QuadraticPenalty.md)
  object.

- beta:

  A numeric vector of length `pen@n_coef`. Unused: no quantity here
  depends on the coefficients.

- theta:

  A named list holding `lambda`. Unused for the same reason.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
a list of one matrix named `lambda`.
[`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md)
a list of one zero matrix named `lambda_lambda`.
[`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md)
a list of one zero vector named `lambda_lambda`.
[`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md)
the single logical `TRUE`.

## Details

[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
returns the stored \\D'PD\\ **in whatever storage the penalty keeps
it**, [`unclass()`](https://rdrr.io/r/base/class.html)ed to strip the
attributes a base matrix may carry and leaving an S4 matrix untouched. A
penalty built with `blocks > 1` therefore answers with a `dgCMatrix`,
and a consumer that needs a base matrix coerces where the two meet.

The zeros are exact. The Hessian is \\\lambda \cdot\\ a constant, so its
second derivative in \\\lambda\\ vanishes identically, and so does the
third derivative of the value in \\\beta\\ and two hyperparameters.

## See also

[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
and its siblings for the generics,
[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
for the branch,
[`zero_pairs()`](https://statmodels7.github.io/penalties7/reference/zero_pairs.md)
for the zeros.

## Examples

``` r
pen <- quadratic_penalty(crossprod(diff(diag(3))))
b <- c(1, -0.5, 0.3)

# The derivative is D'PD, so it does not move with lambda at all.
max(abs(penalty_dhessian(pen, b, list(lambda = 2))$lambda -
          penalty_dhessian(pen, b, list(lambda = 99))$lambda))
#> [1] 0

# And everything above first order is exactly zero.
max(abs(unlist(penalty_d2hessian(pen, b, list(lambda = 2)))))
#> [1] 0
max(abs(unlist(penalty_dcross(pen, b, list(lambda = 2)))))
#> [1] 0
beta_quadratic(pen, list(lambda = 2))
#> [1] TRUE
```
