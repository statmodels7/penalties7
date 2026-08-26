# The Derivative of the Coefficient Hessian in the Hyperparameters

Returns \\\partial S/\partial\theta_m = \partial^3\rho/\partial\beta^2\\
\partial\theta_m\\, one matrix per hyperparameter, saying how the
penalty's curvature in the coefficients moves as each hyperparameter
moves. This is the third-order quantity a marginal likelihood criterion
needs and that the value, the gradient and the Hessian cannot supply
between them.

## Usage

``` r
penalty_dhessian(pen, beta, theta, ...)
```

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- beta:

  A numeric vector of length `pen@n_coef`. Read only by the separable
  branch, whose parent's derivatives depend on the argument; the
  quadratic, additive and structured branches ignore it.

- theta:

  A named list of hyperparameter values, or a named numeric vector
  carrying the same.

- ...:

  Passed to methods. No shipped method reads it.

## Value

A named list of one matrix per hyperparameter, each of side
`pen@n_coef`, named by `pen@params`. The quadratic branch returns the
stored matrix in whatever storage it keeps, so a penalty built with
`blocks > 1` gives a `dgCMatrix`.

## Why a criterion needs it

Estimating hyperparameters by a marginal likelihood means
differentiating a Laplace approximation, whose determinant is \\\lvert
H + S\rvert\\ with \\S\\ the penalty's Hessian in the coefficients.
Differentiating \\\log\lvert H+S\rvert\\ in \\\theta_m\\ needs
\\\partial S/\partial\theta_m\\, and nothing below third order gives it.

## What each branch answers

|  |  |
|----|----|
| branch | \\\partial S/\partial\theta_m\\ |
| [`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md) | \\D'PD\\, free of both \\\lambda\\ and \\\beta\\ |
| [`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md) | the component \\P_k\\ of each smoothing parameter |
| [`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md) | the matrix parameter's own `param_d1` |
| [`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md) | \\-D'\mathrm{diag}(\partial^3\ell/\partial y^2\partial\theta_m)D\\, from [`distributions7::distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.html) |
| [`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md), [`mcp_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md) | rejects |

Every one of the first four is closed form. The base class rejects,
naming what is missing, so a penalty written later either supplies the
three methods and works in a marginal criterion, or says plainly that it
cannot.

A separable penalty with a **kink** rejects as well, whatever its
parent: the third derivative does not exist there, and the mode a
marginal criterion expands around is exactly where a kinked penalty puts
coefficients.

## See also

[`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
for the quantity differentiated,
[`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md)
and
[`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md)
for the other two a marginal criterion asks for,
[`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md)
for the predicate that says whether a third \\\beta\\-derivative is
needed at all.

## Examples

``` r
# A quadratic penalty's Hessian is lambda D'PD, so the derivative is D'PD
# and does not move with lambda.
pen <- quadratic_penalty(crossprod(diff(diag(3))))
b <- c(1, -0.5, 0.3)
penalty_dhessian(pen, b, list(lambda = 2))$lambda
#>      [,1] [,2] [,3]
#> [1,]    1   -1    0
#> [2,]   -1    2   -1
#> [3,]    0   -1    1
identical(penalty_dhessian(pen, b, list(lambda = 2)),
          penalty_dhessian(pen, b, list(lambda = 99)))
#> [1] TRUE

# An additive penalty answers with its components, one per parameter.
add <- additive_penalty(list(diag(3), diag(c(1, 1, 0))))
names(penalty_dhessian(add, b, list(lambda1 = 2, lambda2 = 0.5)))
#> [1] "lambda1" "lambda2"

# A penalty with a kink has no such quantity and says so.
try(penalty_dhessian(lasso_penalty(n_coef = 3), b, list(lambda = 1)))
#> Error : 'separable [fixed laplace2 [mu=0]]' has a kink, so penalty_dhessian() does not exist there and
#>   its hyperparameters cannot be estimated by a marginal criterion.
```
