# The Derivative of the Coefficient Hessian in the Coefficients

Returns \\\partial S/\partial\beta\\\[v\] = \sum_c v_c\\
\partial^3\rho/\partial\beta\\\partial\beta'\\\partial\beta_c\\, the
coefficient Hessian's derivative along a direction \\v\\, as one matrix
of side `pen@n_coef`. It is zero for every penalty quadratic in the
coefficients.

## Usage

``` r
penalty_dhessian_beta(pen, beta, theta, v, ...)
```

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- beta:

  A numeric vector of length `pen@n_coef`.

- theta:

  A named list of hyperparameter values, or a named numeric vector
  carrying the same.

- v:

  A numeric vector of length `pen@n_coef`, the direction.

- ...:

  Passed to methods. No shipped method reads it.

## Value

A square base matrix of side `pen@n_coef`.

## Why a criterion needs it

A marginal criterion reads \\\log\lvert H + S\rvert\\ at the penalized
mode, and the mode moves with the hyperparameters. Where \\S\\ does not
depend on the coefficients that movement reaches the determinant through
\\H\\ alone. Where it does, as for a heavy-tailed prior on a random
effect, it reaches it through \\S\\ as well, and the gradient needs
\\\mathrm{tr}(M\\\partial S/\partial\beta\[v\])\\ with \\v\\ the mode's
movement. A prediction-error criterion reads the same matrix inside the
trace that counts its degrees of freedom.

The direction is contracted here rather than returning the array of
order three, because every consumer reads the array along a direction
and it has \\q^3\\ entries.

## What each branch answers

|  |  |
|----|----|
| branch | \\\partial S/\partial\beta\[v\]\\ |
| [`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md), [`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md), [`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md) | zero |
| [`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md) whose parent is quadratic in the argument | zero |
| [`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md) with a univariate parent otherwise | \\-D'\mathrm{diag}(\ell^{(yyy)}(D\beta)\odot Dv)\\D\\ |
| [`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md) with a multivariate parent otherwise | rejects |
| a kinked parent, [`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md), [`mcp_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md) | rejects |

The univariate row follows from \\S =
-D'\mathrm{diag}(\ell^{(yy)}(D\beta))D\\: differentiating
\\\ell^{(yy)}((D\beta)\_j)\\ in \\\beta_c\\ gives
\\\ell^{(yyy)}((D\beta)\_j) D\_{jc}\\, and summing against \\v_c\\ gives
\\(Dv)\_j\\. The third response derivative is read from
[`distributions7::distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.html),
which is closed form for every location family and so for a Student t
prior.

A multivariate parent that is not quadratic, a multivariate t prior
among them, would need the third response derivative as an array per
block, which distributions7 does not supply, so it rejects rather than
returning a matrix missing that piece. Whether the parent is quadratic
is asked of
[`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md)
first, so a Gaussian prior of any dimension answers zero without
reaching the parent at all.

## See also

[`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
for the quantity differentiated,
[`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md)
for the predicate that says it is zero,
[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
for the derivative in the hyperparameters.

## Examples

``` r
b <- c(1, -0.5, 0.3)
v <- c(0.2, 0.1, -0.4)

# Zero for a quadratic penalty, whose Hessian does not move with beta.
penalty_dhessian_beta(quadratic_penalty(diag(3)), b, list(lambda = 2), v)
#>      [,1] [,2] [,3]
#> [1,]    0    0    0
#> [2,]    0    0    0
#> [3,]    0    0    0

# A Student t prior's Hessian does move, and this is its derivative along v.
h <- heavy_penalty(n_coef = 3)
th <- list(sigma = 1, nu = 4)
D <- penalty_dhessian_beta(h, b, th, v)
eps <- 1e-6
num <- (penalty_hessian(h, b + eps * v, th) -
        penalty_hessian(h, b - eps * v, th)) / (2 * eps)
max(abs(D - num))
#> [1] 8.447346e-11

# A kinked parent has no such derivative at the kink and says so.
try(penalty_dhessian_beta(lasso_penalty(n_coef = 3), b, list(lambda = 1), v))
#> Error : 'separable [fixed laplace2 [mu=0]]' has a kink, so penalty_dhessian_beta() does not exist there and
#>   its hyperparameters cannot be estimated by a marginal criterion.
```
