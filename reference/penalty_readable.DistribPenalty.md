# Interpretable Quantities of a Separable Penalty

Returns `NULL` for a univariate parent, whose hyperparameters are
already interpretable quantities, and the result of
[`distributions7::mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.html)
for a multivariate one, whose hyperparameters are the free values of a
matrix parameter.

## Arguments

- penalty:

  A
  [`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
  object.

- theta:

  A named list of the parent's free parameters.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

`NULL` when `penalty@block` is `1L`. Otherwise the method returns the
list that
[`distributions7::mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.html)
produces, with the elements `value`, `jacobian`, `transform` and
`block`.

## Details

The hyperparameters of a univariate parent are quantities such as a
scale, a rate or a shape, each on its own scale, so there is nothing to
derive and the method returns `NULL`. Those of a multivariate parent are
coordinates such as log-Cholesky entries, which are not directly
interpretable; the interpretable quantities are the standard deviations
and the correlations of the effects within a block, and the parent
provides them through
[`distributions7::mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.html).

## See also

[`penalty_readable()`](https://statmodels7.github.io/penalties7/reference/penalty_readable.md)
for the generic,
[`distributions7::mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.html)
for the quantities that a multivariate family provides.

## Examples

``` r
# For a univariate parent the rate is already the quantity of interest, so
# the result is NULL.
penalty_readable(lasso_penalty(n_coef = 3), list(lambda = 1))
#> NULL

# A bivariate Gaussian prior has two standard deviations and a correlation.
mv <- distrib_penalty(
  distributions7::fixed(distributions7::mvgaussian1_distrib(2),
                        mu1 = 0, mu2 = 0), n_coef = 6)
r <- penalty_readable(mv, list(sigma_log_L1 = 0.2, sigma_log_L2 = -0.1,
                               sigma_L2.1 = 0.5))
r$value
#>     sd_v1     sd_v2 cor_v1_v2 
#> 1.2214028 1.0337943 0.4836552 
r$transform
#>     sd_v1     sd_v2 cor_v1_v2 
#>     "log"     "log"   "atanh" 
```
