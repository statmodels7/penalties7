# The Derivative of the Coefficient Hessian in the Hyperparameters

Returns \\\partial S/\partial\theta_m = \partial^3\rho/\partial\beta^2\\
\partial\theta_m\\, one matrix per hyperparameter, which describes how
the curvature of the penalty in the coefficients changes with each
hyperparameter. A marginal likelihood criterion needs this third-order
quantity, and the value, the gradient and the Hessian do not determine
it.

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

  A numeric vector of length `pen@n_coef`; the generic signals an error
  for any other length. Only the separable branch reads its values,
  because the derivatives of the parent depend on them; the quadratic,
  additive and structured branches ignore them.

- theta:

  A named list of hyperparameter values, or a named numeric vector
  carrying the same. The generic checks that every hyperparameter is
  present as a single number inside its bounds.

- ...:

  Passed to methods. The shipped methods do not read it.

## Value

A named list of one matrix per hyperparameter, each of side
`pen@n_coef`, named by `pen@params`. The quadratic branch returns the
stored matrix in the storage that the penalty keeps, so a penalty built
with `blocks > 1` gives a `dgCMatrix`; the other branches return base
matrices.

## Role in a marginal criterion

Estimating hyperparameters by a marginal likelihood means
differentiating a Laplace approximation, whose determinant is \\\lvert
H + S\rvert\\ with \\S\\ the Hessian of the penalty in the coefficients.
The derivative of \\\log\lvert H+S\rvert\\ in \\\theta_m\\ contains
\\\partial S/\partial\theta_m\\, which is a third derivative of
\\\rho\\. The derivatives are taken in \\\theta\\ on its own scale,
without a link.

## The value on each branch

|  |  |
|----|----|
| branch | \\\partial S/\partial\theta_m\\ |
| [`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md) | \\D'PD\\, free of both \\\lambda\\ and \\\beta\\ |
| [`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md) | the component \\D'P_kD\\ of each smoothing parameter, which is \\P_k\\ under the identity map |
| [`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md) | the derivative `param_d1` of the matrix parameter |
| [`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md), univariate parent | \\-D'\mathrm{diag}(\partial^3\ell/\partial y^2\partial\theta_m)D\\, from [`distributions7::distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.html) |
| [`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md), multivariate parent | the same with a block-diagonal middle matrix |
| [`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md), [`mcp_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md) | an error, from the base class |

The quadratic and additive rows are exact. The structured row is exact
for the structures that parameters7 ships, and the separable rows are
exact when the parent supplies
[`distributions7::distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.html)
in closed form; otherwise the numerical fallbacks of those packages
apply. On the base class the method signals an error that names the
penalty, so a penalty class written later is usable in a marginal
criterion only if it supplies `penalty_dhessian()`,
[`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md)
and
[`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md).

A separable penalty with a kink signals an error as well, whatever its
parent. The third derivative does not exist at the kink, and a kinked
penalty places coefficients at the mode at which a marginal criterion is
expanded.

## See also

[`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
for the quantity differentiated,
[`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md)
and
[`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md)
for the other two quantities that a marginal criterion requires,
[`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md)
for the predicate that indicates whether a third \\\beta\\-derivative is
needed.

## Examples

``` r
# A quadratic penalty's Hessian is lambda D'PD, so the derivative is D'PD
# and does not change with lambda.
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

# An additive penalty returns its components, one per parameter.
add <- additive_penalty(list(diag(3), diag(c(1, 1, 0))))
names(penalty_dhessian(add, b, list(lambda1 = 2, lambda2 = 0.5)))
#> [1] "lambda1" "lambda2"

# A penalty with a kink has no such quantity, and the call is an error.
try(penalty_dhessian(lasso_penalty(n_coef = 3), b, list(lambda = 1)))
#> Error : 'separable [fixed laplace2 [mu=0]]' has a kink, so penalty_dhessian() does not exist there and
#>   its hyperparameters cannot be estimated by a marginal criterion.
```
