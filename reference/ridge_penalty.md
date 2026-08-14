# Named Separable Penalties

The canonical instances of
[`distrib_penalty`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md),
shipped as constructors so the model layer can name what it means. Each
is written on the chart whose hyperparameter MEASURES THE SHRINKAGE, so
that a larger value shrinks harder in all of them. Ridge is the
exception to the branch rather than to the rule: it is the Gaussian
prior at zero, and that prior written by its PRECISION is exactly the
quadratic penalty at the identity, the same value to the last bit, so it
is built there and its hyperparameter is the lambda that branch already
carries – one name for one number. The lasso is the Laplace in location
and rate (\`laplace2\`) at zero, so the free hyperparameter is the rate
\\\lambda\\ and the value is \\\lambda\lVert D\beta\rVert_1\\ up to its
constant, with the kink declared; the heavy-tailed prior is the Student
t at zero, whose \\\nu\\ is estimable exactly because the normalizing
constant is kept. The elastic net is the product of the Laplace and the
Gaussian at zero, normalized
([`enet_distrib`](https://statmodels7.github.io/distributions7/reference/enet_distrib.html)),
so its hyperparameters are the overall rate \\\lambda\\ and the mixing
weight \\\alpha\\ and its value is \\\lambda\\\alpha\lVert
D\beta\rVert_1 + (1-\alpha)\lVert D\beta\rVert_2^2/2\\\\ up to a
constant. That constant depends on both hyperparameters, which is what
makes them estimable by a marginal criterion and what a penalty written
as a formula would not have.

## Usage

``` r
ridge_penalty(map = NULL, n_coef = 1L)

lasso_penalty(map = NULL, n_coef = 1L)

elasticnet_penalty(map = NULL, n_coef = 1L)

heavy_penalty(map = NULL, n_coef = 1L)
```

## Arguments

- map:

  The matrix \\D\\, or `NULL` (default) for the identity.

- n_coef:

  The number of coefficients when `map` is `NULL`.

## Value

An object of class `DistribPenalty`.

## References

Hoerl, A. E. and Kennard, R. W. (1970). Ridge regression: biased
estimation for nonorthogonal problems. *Technometrics* 12, 55-67.

Tibshirani, R. (1996). Regression shrinkage and selection via the lasso.
*Journal of the Royal Statistical Society, Series B* 58, 267-288.

Zou, H. and Hastie, T. (2005). Regularization and variable selection via
the elastic net. *Journal of the Royal Statistical Society, Series B*
67, 301-320.

## See also

[`distrib_penalty`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md),
[`scad_penalty`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md),
[`quadratic_penalty`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)

## Examples

``` r
pen <- ridge_penalty(n_coef = 2)
penalty_gradient(pen, c(1, -1), list(lambda = 1))
#> [1]  1 -1
```
