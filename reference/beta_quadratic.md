# Is a Penalty Quadratic in the Coefficients?

`TRUE` when \\\partial^3\rho/\partial\beta^3\\ is exactly zero, so a
consumer can skip that third derivative altogether. Asked of the penalty
rather than measured by a consumer, which would be guessing at a
property the penalty knows.

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
  carrying the same. Read by the separable branch, to evaluate the
  parent; the other branches answer a constant.

- ...:

  Passed to methods. No shipped method reads it.

## Value

A single logical. `FALSE` from the base class, so a penalty that says
nothing is treated as needing the third derivative.

## Three independent properties

This is not
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md),
which is about the whole construction. The three questions a consumer
may ask are independent:

|  |  |  |  |
|----|----|----|----|
| penalty | [`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md) | `beta_quadratic()` | \\S\\ linear in \\\theta\\ |
| [`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md) | `TRUE` | `TRUE` | yes |
| [`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md) | `TRUE` | `TRUE` | no |
| [`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md) | `FALSE` | `TRUE` | yes |
| a Gaussian [`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md) | `FALSE` | `TRUE` | no |
| [`heavy_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md) | `FALSE` | `FALSE` | no |

## How the separable branch answers

By probing, not by differentiating three times. The log-density is
quadratic in the response exactly when its second derivative there does
not depend on it, so
[`distributions7::distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.html)
is read at four points and compared for constancy at a tolerance of
`1e-12`. The question is put to the second derivative because that one
is analytic for almost every family, where the third is often a finite
difference whose noise no threshold separates from a true zero.

## What `TRUE` does not promise

It does not say the marginal derivatives are available. A lasso answers
`TRUE`, its log-density being linear in the response away from the kink,
and
[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
still rejects for it: the third derivative is zero where it exists and
undefined at the kink, which is where a selecting fit puts its
coefficients.

## See also

[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
for the different question about the construction,
[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
for what a consumer skips when this is `TRUE`,
[`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
for the derivative that is then constant.

## Examples

``` r
# The two branches that answer TRUE to both predicates.
beta_quadratic(quadratic_penalty(diag(3)), list(lambda = 1))
#> [1] TRUE
beta_quadratic(structured_penalty(
  parameters7::log_cholesky(2, role = "precision")),
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
