# Second Derivatives of the Coefficient Hessian in the Coefficients and the Hyperparameters

These two functions return the derivatives of
[`penalty_dhessian_beta()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian_beta.md)
that the second derivative of a marginal criterion contains.
`penalty_d2hessian_beta()` is \\\partial^2 S/\partial\beta^2\\\[v,
w\]\\, the second derivative of the Hessian in the coefficients
contracted along two directions. `penalty_dhessian_beta_theta()` is
\\\partial^2 S/\partial\beta\\\partial\theta_m\\\[v\]\\, the derivative
of \\\partial S/\partial\beta\\\[v\]\\ in each hyperparameter.

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

  A numeric vector of length `pen@n_coef`; the generics signal an error
  for any other length.

- theta:

  A named list of hyperparameter values, or a named numeric vector
  carrying the same, checked by the generics as in
  [`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md).

- v, w:

  Numeric vectors of length `pen@n_coef`, the directions; the generics
  signal an error for any other length.

- ...:

  Passed to methods. The shipped methods do not read it.

## Value

`penalty_d2hessian_beta()` returns a square base matrix of side
`pen@n_coef`; `penalty_dhessian_beta_theta()` a list of such matrices,
one per hyperparameter, keyed by `pen@params`.

## Details

Write \\S_m = \partial S/\partial\theta_m\\, \\T\[v\] = \partial
S/\partial\beta\\\[v\]\\, and \\b_m\\ and \\b\_{ml}\\ for the first and
second derivatives of the mode in the hyperparameters. The matrix \\K_m
= S_m + T\[b_m\]\\ enters the gradient of the criterion, and its
derivative in \\\theta_l\\ has three parts that are present only for a
penalty whose Hessian depends on the coefficients: \\T\[b\_{ml}\]\\,
from
[`penalty_dhessian_beta()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian_beta.md);
the second derivative in the coefficients along \\b_m\\ and \\b_l\\,
from `penalty_d2hessian_beta()`; and the two matrices \\\partial
S_m/\partial\beta\\\[b_l\]\\ and \\\partial
S_l/\partial\beta\\\[b_m\]\\. Both of these come from
`penalty_dhessian_beta_theta()`, the first as its component \\m\\ along
\\b_l\\ and the second as its component \\l\\ along \\b_m\\, because by
the symmetry of mixed partials in \\\beta\\ and \\\theta\\ the
derivative of \\T\[b_m\]\\ in \\\theta_l\\ equals \\\partial
S_l/\partial\beta\\\[b_m\]\\. The two matrices are in general different.

## The value on each branch

|  |  |  |
|----|----|----|
| branch | \\\partial^2 S/\partial\beta^2\[v,w\]\\ | \\\partial^2 S/\partial\beta\partial\theta_m\[v\]\\ |
| quadratic, additive, structured | zero | zero for every hyperparameter |
| separable, parent quadratic in its argument | zero | zero for every hyperparameter |
| separable, univariate parent otherwise | \\-D'\mathrm{diag}(\ell^{(yyyy)} \odot Dv \odot Dw)D\\ | \\-D'\mathrm{diag}(\partial\_{\theta_m}\ell^{(yyy)} \odot Dv)D\\ |
| separable, multivariate parent otherwise | the same per block, contracted by [`dp_contract()`](https://statmodels7.github.io/penalties7/reference/dp_contract.md) | the same per block |
| a kinked parent, [`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md), [`mcp_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md) | an error | an error |

The univariate rows follow from \\\partial S/\partial\beta\[v\] =
-D'\mathrm{diag}(\ell^{(yyy)}(D\beta) \odot Dv)D\\ by differentiating
\\\ell^{(yyy)}((D\beta)\_j)\\ once more, in \\\beta\\ along \\w\\ or in
\\\theta_m\\. The fourth response derivative is
[`distributions7::distrib_deriv4_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.html),
closed form for every family that distributions7 ships and one stencil
for a family defined by its density alone. The mixed one is
[`distributions7::distrib_cross3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross3_y.html),
closed form for the location families listed on its page, the Student t
among them, and one central difference of the third response derivative
for the other families. A multivariate parent supplies the two as arrays
per block, of \\p^4\\ and \\p^3\\ entries, which
[`dp_contract()`](https://statmodels7.github.io/penalties7/reference/dp_contract.md)
contracts against the coordinates of the block in the directions; the
multivariate Student t supplies both. The derivatives are taken in
\\\theta\\ on its own scale.

The methods of the base class and of each branch are described on
[`penalty_d2hessian_beta.penalty()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.penalty.md).

## See also

[`penalty_dhessian_beta()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian_beta.md)
for the quantity differentiated,
[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
for the derivative of the Hessian in the hyperparameters,
[`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md)
for the predicate under which both are zero.

## Examples

``` r
h <- heavy_penalty(n_coef = 3)
b <- c(1, -0.5, 0.3)
v <- c(0.2, 0.1, -0.4)
w <- c(-0.3, 0.5, 0.1)
th <- list(sigma = 1, nu = 4)

# The second derivative along two directions is compared with a difference
# of the first derivative along one of them.
eps <- 1e-6
D2 <- penalty_d2hessian_beta(h, b, th, v, w)
num <- (penalty_dhessian_beta(h, b + eps * w, th, v) -
        penalty_dhessian_beta(h, b - eps * w, th, v)) / (2 * eps)
max(abs(D2 - num))
#> [1] 2.504424e-11

# The movement in the degrees of freedom is compared with a central
# difference.
Dnu <- penalty_dhessian_beta_theta(h, b, th, v)$nu
num <- (penalty_dhessian_beta(h, b, list(sigma = 1, nu = 4 + eps), v) -
        penalty_dhessian_beta(h, b, list(sigma = 1, nu = 4 - eps), v)) /
  (2 * eps)
max(abs(Dnu - num))
#> [1] 7.828853e-11

# The result is zero for a quadratic penalty.
penalty_d2hessian_beta(quadratic_penalty(diag(3)), b, list(lambda = 2), v, w)
#>      [,1] [,2] [,3]
#> [1,]    0    0    0
#> [2,]    0    0    0
#> [3,]    0    0    0
```
