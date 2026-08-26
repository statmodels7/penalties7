# The Named Penalties

The four penalties a modeling layer names directly: `ridge_penalty()`,
`lasso_penalty()`, `elasticnet_penalty()` and `heavy_penalty()`. Each is
a particular prior centered at zero, and each is written on the chart
whose hyperparameter **measures the shrinkage**, so that a larger value
shrinks harder in all four.

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
  default) for the identity. Given a map, the number of coefficients is
  `ncol(map)` and `n_coef` is ignored.

- n_coef:

  The number of coefficients when `map` is `NULL`. A single whole
  number, `1L` by default.

## Value

`ridge_penalty()` a
[`QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/QuadraticPenalty.md)
with the hyperparameter `lambda`, a precision, on \\(0, \infty)\\.
`lasso_penalty()` a
[`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
with `lambda`, a rate, on \\(0, \infty)\\ and a kink at zero.
`elasticnet_penalty()` a
[`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
with `lambda` on \\(0, \infty)\\ and `alpha` on \\(0, 1)\\, and a kink
at zero. `heavy_penalty()` a
[`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
with `sigma` and `nu` both on \\(0, \infty)\\, and no kink.

## Ridge, which is on the other branch

`ridge_penalty()` is the Gaussian prior at zero, and that prior written
by its **precision** is exactly
[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
at the identity matrix, the same value to the last bit. It is therefore
built there, and its hyperparameter is the `lambda` that branch already
carries: one name for one number. It is the only one of the four that is
not a `DistribPenalty`, and
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
answers `TRUE` for it alone.

The separable twin still exists and is one line,
`distrib_penalty(fixed(gaussian1_distrib(), mu = 0), n_coef = q)`, whose
hyperparameter is the standard deviation \\\sigma\\ instead. Reach for
it when a standard deviation is the quantity you want reported.

## Lasso

The Laplace in location and rate (`laplace2`) held at zero, so the free
hyperparameter is the rate \\\lambda\\ and

\$\$\rho(\beta;\lambda) = \lambda\lVert D\beta\rVert_1 -
q\log(\lambda/2),\$\$

with \\q\\ the number of values the penalty is read at. The kink at zero
is declared, so
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
reports it and
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
gives the soft threshold.

## Elastic net

The product of the Laplace and the Gaussian at zero, normalized
([`distributions7::enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.html)),
so the hyperparameters are the overall rate \\\lambda\\ and the mixing
weight \\\alpha\\ on \\(0, 1)\\:

\$\$\rho(\beta;\lambda,\alpha) = \lambda\left\\ \alpha\lVert
D\beta\rVert_1 + (1-\alpha)\lVert D\beta\rVert_2^2/2\right\\ +
c(\lambda, \alpha).\$\$

The constant depends on **both** hyperparameters, which is why a
marginal criterion can estimate them: measured at \\(\lambda, \alpha)\\
of \\(1, 0.5)\\, \\(2, 0.5)\\ and \\(1, 0.9)\\, the per-coordinate
constant is `0.7805`, `0.2711` and `0.7001`. A penalty written as a
formula, with no constant, would have neither hyperparameter identified.

## The heavy-tailed prior

The Student t at zero, with a free scale \\\sigma\\ and a free \\\nu\\.
Its \\\nu\\ is estimable exactly because the normalizing constant is
kept, and estimating it is the point: a small \\\nu\\ shrinks small
coefficients like a Gaussian prior and leaves large ones nearly alone,
which a Gaussian cannot do at any scale. Alone among the four it
declares no kink, so it produces no exact zeros.

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
for the construction behind three of them,
[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
for the one behind ridge,
[`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
and
[`mcp_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
for the non-convex alternatives,
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
for the operators these have.

## Examples

``` r
b <- c(1, -0.5, 0.3)

# Ridge is on the quadratic branch, so its lambda is a precision.
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

# Only the lasso and the elastic net have a kink, so only they can set a
# coefficient exactly to zero.
penalty_kinks(l, list(lambda = 1))
#> [1] 0
penalty_kinks(heavy_penalty(n_coef = 3), list(sigma = 1, nu = 4))
#> numeric(0)

# The heavy-tailed prior shrinks a large coefficient far less than a
# ridge at a comparable scale, and the small ones about as much.
penalty_prox(heavy_penalty(n_coef = 3), c(4, 0.3, -1), 1,
             list(sigma = 1, nu = 3))
#> [1]  3.0000000  0.1289779 -0.4442398
penalty_prox(r, c(4, 0.3, -1), 1, list(lambda = 1))
#> [1]  2.00  0.15 -0.50
```
