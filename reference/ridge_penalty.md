# The Named Penalties

These functions construct the four penalties that a modeling layer names
directly: `ridge_penalty()`, `lasso_penalty()`, `elasticnet_penalty()`
and `heavy_penalty()`. Each is a prior centered at zero. In the first
three the hyperparameter `lambda` is a precision or a rate, and a larger
value shrinks harder. In `heavy_penalty()` the hyperparameter `sigma` is
a scale, and a larger value shrinks less.

## Usage

``` r
ridge_penalty(map = NULL, n_coef = 1L)

lasso_penalty(map = NULL, n_coef = 1L)

elasticnet_penalty(map = NULL, n_coef = 1L)

heavy_penalty(map = NULL, n_coef = 1L)
```

## Arguments

- map:

  The matrix \\D\\, with one column per coefficient, or `NULL` (the
  default) for the identity. If a map is given, the number of
  coefficients is `ncol(map)` and `n_coef` is ignored.

- n_coef:

  The number of coefficients when `map` is `NULL`. It is a single whole
  number, and the default is `1L`.

## Value

`ridge_penalty()` returns a
[`QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/QuadraticPenalty.md)
with the hyperparameter `lambda`, a precision, on \\(0, \infty)\\.
`lasso_penalty()` returns a
[`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
with `lambda`, a rate, on \\(0, \infty)\\ and a kink at zero.
`elasticnet_penalty()` returns a
[`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
with `lambda` on \\(0, \infty)\\ and `alpha` on \\(0, 1)\\, and a kink
at zero. `heavy_penalty()` returns a
[`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
with `sigma` and `nu` both on \\(0, \infty)\\, and no kink.

## Ridge

`ridge_penalty()` is the Gaussian prior at zero written by its
precision, which is
[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
at the identity matrix: its value agrees with that of the separable
Gaussian prior at \\\sigma = \lambda^{-1/2}\\ up to rounding. It is
therefore built on the quadratic branch, and its hyperparameter is the
`lambda` of that branch. It is the only one of the four that is not a
`DistribPenalty`, and the only one for which
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
returns `TRUE`.

The separable form is
`distrib_penalty(fixed(gaussian1_distrib(), mu = 0), n_coef = q)`, whose
hyperparameter is the standard deviation \\\sigma\\ instead. It is the
form to use when a standard deviation is the quantity to be reported.

## Lasso

The lasso is the Laplace in location and rate (`laplace2`) held at zero,
so the free hyperparameter is the rate \\\lambda\\ and

\$\$\rho(\beta;\lambda) = \lambda\lVert D\beta\rVert_1 -
m\log(\lambda/2),\$\$

where \\m\\ is the number of values at which the penalty is read (the
number of rows of \\D\\, or the number of coefficients without a map).
Under a square nonsingular map the value also carries \\-\log\lvert\det
D\rvert\\, as described in
[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md).
The kink at zero is declared, so
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
reports it and
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
gives the soft threshold.

## Elastic net

The elastic net is the product of the Laplace and the Gaussian at zero,
normalized
([`distributions7::enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.html)),
so the hyperparameters are the overall rate \\\lambda\\ and the mixing
weight \\\alpha\\ on \\(0, 1)\\:

\$\$\rho(\beta;\lambda,\alpha) = \lambda\left\\ \alpha\lVert
D\beta\rVert_1 + (1-\alpha)\lVert D\beta\rVert_2^2/2\right\\ +
m\\c(\lambda, \alpha),\$\$

where \\m\\ is as above and \\c(\lambda,\alpha) = \log\int
\exp\[-\lambda\\\alpha\lvert x\rvert + (1-\alpha)x^2/2\\\]\\dx\\ is the
constant of one coordinate. The same Jacobian term is added under a
square nonsingular map. The constant depends on both hyperparameters, so
a marginal criterion can estimate them. Without the constant, both
hyperparameters would be unidentified.

## Heavy-tailed prior

The heavy-tailed prior is the Student t at zero, with a free scale
\\\sigma\\ and free degrees of freedom \\\nu\\. The normalizing constant
is kept, which makes \\\nu\\ estimable by a marginal criterion. With a
small \\\nu\\ the prior shrinks small coefficients as a Gaussian prior
does and leaves large ones nearly unchanged, which a Gaussian prior does
not do at any scale. Like the ridge, it has no kink, so it does not
produce exact zeros.

## References

Hoerl, A. E. and Kennard, R. W. (1970). Ridge regression: biased
estimation for nonorthogonal problems. *Technometrics* **12**, 55-67.

Tibshirani, R. (1996). Regression shrinkage and selection via the lasso.
*Journal of the Royal Statistical Society, Series B* **58**, 267-288.

Zou, H. and Hastie, T. (2005). Regularization and variable selection via
the elastic net. *Journal of the Royal Statistical Society, Series B*
**67**, 301-320.

## See also

[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md)
for the construction of the lasso, the elastic net and the heavy-tailed
prior,
[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
for the construction of the ridge,
[`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
and
[`mcp_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
for the non-convex alternatives,
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
for their proximal operators.

## Examples

``` r
b <- c(1, -0.5, 0.3)

# The ridge is on the quadratic branch, and its lambda is a precision.
r <- ridge_penalty(n_coef = 3)
class(r)[1]
#> [1] "penalties7::QuadraticPenalty"
is_quadratic(r)
#> [1] TRUE
penalty_gradient(r, b, list(lambda = 2))
#> [1]  2.0 -1.0  0.6

# The lasso's value is lambda times the L1 norm, plus its constant.
l <- lasso_penalty(n_coef = 3)
penalty_value(l, b, list(lambda = 0.5))
#> [1] 5.058883
0.5 * sum(abs(b)) - 3 * log(0.5 / 2)
#> [1] 5.058883

# Of these four, only the lasso and the elastic net have a kink, so only
# they can set a coefficient exactly to zero.
penalty_kinks(l, list(lambda = 1))
#> [1] 0
penalty_kinks(heavy_penalty(n_coef = 3), list(sigma = 1, nu = 4))
#> numeric(0)

# The heavy-tailed prior shrinks the large coefficient far less than a
# unit ridge, and the small ones somewhat more.
penalty_prox(heavy_penalty(n_coef = 3), c(4, 0.3, -1), 1,
             list(sigma = 1, nu = 3))
#> [1]  3.0000000  0.1289779 -0.4442398
penalty_prox(r, c(4, 0.3, -1), 1, list(lambda = 1))
#> [1]  2.00  0.15 -0.50
```
