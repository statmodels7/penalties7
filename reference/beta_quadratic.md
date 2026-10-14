# Is a Penalty Quadratic in the Coefficients?

`TRUE` when \\\partial^3\rho/\partial\beta^3\\ is exactly zero, so that
a consumer can skip that third derivative. The value comes from a method
of the class of the penalty, so a consumer does not have to test the
property numerically.

## Usage

``` r
beta_quadratic(pen, theta, ...)
```

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- theta:

  A named list of hyperparameter values, or a named numeric vector
  carrying the same, holding every name in `pen@params`. The generic
  checks it as
  [`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
  does, so a missing hyperparameter or a value outside its open bounds
  is an error on every branch. The separable branch evaluates the parent
  at it, and the other branches return a constant without reading its
  values.

- ...:

  Passed to methods. No shipped method reads it.

## Value

A single logical. The base class returns `FALSE`, so a penalty class
without a method of its own is treated as needing the third derivative.

## Three independent properties

This predicate differs from
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md),
which reports whether the penalty is a quadratic form whose matrix, rank
and log pseudo-determinant are available through
[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
and its siblings. The two predicates and the linearity of \\S\\ in
\\\theta\\ are independent properties:

|  |  |  |  |
|----|----|----|----|
| penalty | [`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md) | `beta_quadratic()` | \\S\\ linear in \\\theta\\ |
| [`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md) | `TRUE` | `TRUE` | yes |
| [`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md) | `TRUE` | `TRUE` | no |
| [`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md) | `TRUE` | `TRUE` | yes |
| a Gaussian [`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md) | `FALSE` | `TRUE` | no |
| [`heavy_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md) | `FALSE` | `FALSE` | no |

## The method of the separable branch

The log-density is quadratic in the response exactly when its second
derivative in the response does not depend on the response. The method
therefore evaluates
[`distributions7::distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.html)
of the parent at four points and tests whether the values agree across
the points with [`all.equal()`](https://rdrr.io/r/base/all.equal.html)
at a tolerance of `1e-12`, which is a relative tolerance. No third
derivative is computed. The page
[`beta_quadratic.DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.DistribPenalty.md)
describes the points and the multivariate case.

## The meaning of `TRUE`

`TRUE` does not imply that the marginal derivatives are available. A
lasso gives `TRUE`, because its log-density is linear in the response
away from the kink, and
[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
still signals an error for it: the third derivative is zero where it
exists and undefined at the kink, which is where a selecting fit places
its coefficients.

## See also

[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
for the property of the construction,
[`penalty_dhessian_beta()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian_beta.md)
for the derivative that is zero when this is `TRUE`,
[`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
for the Hessian that is then constant in the coefficients.

## Examples

``` r
# Two of the three branches that give TRUE for both predicates.
beta_quadratic(quadratic_penalty(diag(3)), list(lambda = 1))
#> [1] TRUE
beta_quadratic(structured_penalty(
  parameters7::log_cholesky(2)),
  list(log_L1 = 0, log_L2 = 0, L2.1 = 0))
#> [1] TRUE

# A Gaussian prior is quadratic in beta and is not a quadratic penalty.
g <- distrib_penalty(
  distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0),
  n_coef = 3)
c(beta_quadratic(g, list(sigma = 1)), is_quadratic(g))
#> [1]  TRUE FALSE

# A Student t prior is neither.
beta_quadratic(heavy_penalty(n_coef = 3), list(sigma = 1, nu = 4))
#> [1] FALSE

# And TRUE does not mean the marginal derivatives exist.
beta_quadratic(lasso_penalty(n_coef = 3), list(lambda = 1))
#> [1] TRUE
try(penalty_dhessian(lasso_penalty(n_coef = 3), c(1, 0, -1),
                     list(lambda = 1)))
#> Error : 'separable [fixed laplace2 [mu=0]]' has a kink, so penalty_dhessian() does not exist there and
#>   its hyperparameters cannot be estimated by a marginal criterion.
```
