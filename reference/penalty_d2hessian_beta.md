# How the Movement of the Coefficient Hessian Moves

The two derivatives of
[`penalty_dhessian_beta()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian_beta.md)
that the second derivative of a marginal criterion reads.
`penalty_d2hessian_beta()` is \\\partial^2 S/\partial\beta^2\\\[v,
w\]\\, the Hessian's second derivative in the coefficients contracted
along two directions. `penalty_dhessian_beta_theta()` is \\\partial^2
S/\partial\beta\\\partial\theta_m\\\[v\]\\, the derivative of \\\partial
S/\partial\beta\\\[v\]\\ in each hyperparameter.

## Usage

``` r
penalty_d2hessian_beta(pen, beta, theta, v, w, ...)

penalty_dhessian_beta_theta(pen, beta, theta, v, ...)
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

- v, w:

  Numeric vectors of length `pen@n_coef`, the directions.

- ...:

  Passed to methods. No shipped method reads it.

## Value

`penalty_d2hessian_beta()` a square base matrix of side `pen@n_coef`;
`penalty_dhessian_beta_theta()` a list of such matrices, one per
hyperparameter, keyed by `pen@params`.

## Details

Differentiating the criterion's gradient once more in the
hyperparameters, the determinant's matrix \\K_m = S_m + T\[b_m\]\\ moves
in three ways that a penalty whose Hessian depends on the coefficients
adds to: through \\T\[b\_{ml}\]\\, which
[`penalty_dhessian_beta()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian_beta.md)
supplies; through the second derivative in the coefficients along the
two directions the mode moves in; and through \\\partial
S_m/\partial\beta\\\[b_l\]\\ and \\\partial
S_l/\partial\beta\\\[b_m\]\\, which are one quantity by the symmetry of
mixed partials.

## What each branch answers

|  |  |  |
|----|----|----|
| branch | \\\partial^2 S/\partial\beta^2\[v,w\]\\ | \\\partial^2 S/\partial\beta\partial\theta_m\[v\]\\ |
| quadratic, additive, structured | zero | zero for every hyperparameter |
| separable, parent quadratic in its argument | zero | zero for every hyperparameter |
| separable, univariate parent otherwise | \\-D'\mathrm{diag}(\ell^{(yyyy)} \odot Dv \odot Dw)D\\ | \\-D'\mathrm{diag}(\partial\_{\theta_m}\ell^{(yyy)} \odot Dv)D\\ |
| separable, multivariate parent otherwise | rejects | rejects |
| a kinked parent, [`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md), [`mcp_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md) | rejects | rejects |

The univariate rows follow from \\\partial S/\partial\beta\[v\] =
-D'\mathrm{diag}(\ell^{(yyy)}(D\beta) \odot Dv)D\\ by differentiating
\\\ell^{(yyy)}((D\beta)\_j)\\ once more, in \\\beta\\ along \\w\\ or in
\\\theta_m\\. The fourth response derivative is
[`distributions7::distrib_deriv4_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.html)
and the mixed one
[`distributions7::distrib_cross3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross3_y.html),
both closed for every location family and so for a Student t prior.

## See also

[`penalty_dhessian_beta()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian_beta.md)
for the quantity differentiated,
[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
for the derivative of the Hessian in the hyperparameters,
[`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md)
for the predicate that says both are zero.

## Examples

``` r
h <- heavy_penalty(n_coef = 3)
b <- c(1, -0.5, 0.3)
v <- c(0.2, 0.1, -0.4)
w <- c(-0.3, 0.5, 0.1)
th <- list(sigma = 1, nu = 4)

# The second derivative along two directions, against a difference of the
# first along one of them.
eps <- 1e-6
D2 <- penalty_d2hessian_beta(h, b, th, v, w)
num <- (penalty_dhessian_beta(h, b + eps * w, th, v) -
        penalty_dhessian_beta(h, b - eps * w, th, v)) / (2 * eps)
max(abs(D2 - num))
#> [1] 1.363723e-11

# The movement in the degrees of freedom.
Dnu <- penalty_dhessian_beta_theta(h, b, th, v)$nu
num <- (penalty_dhessian_beta(h, b, list(sigma = 1, nu = 4 + eps), v) -
        penalty_dhessian_beta(h, b, list(sigma = 1, nu = 4 - eps), v)) /
  (2 * eps)
max(abs(Dnu - num))
#> [1] 7.828854e-11

# Zero for a quadratic penalty.
penalty_d2hessian_beta(quadratic_penalty(diag(3)), b, list(lambda = 2), v, w)
#>      [,1] [,2] [,3]
#> [1,]    0    0    0
#> [2,]    0    0    0
#> [3,]    0    0    0
```
