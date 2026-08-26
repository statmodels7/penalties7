# What a Separable Penalty's Hyperparameters Are About

Returns `NULL` for a univariate parent, whose hyperparameters are
already the quantities a reader reads, and
[`distributions7::mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.html)
for a multivariate one, whose hyperparameters are the free values of a
matrix parameter.

## Arguments

- pen:

  A
  [`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
  object.

- theta:

  A named list of the parent's free parameters.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

`NULL` when `pen@block` is `1L`. Otherwise the list
[`distributions7::mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.html)
returns: `value`, `jacobian`, `transform` and `block`.

## Details

A univariate parent's hyperparameters are a scale, a rate, a shape, each
read on its own scale already, so there is nothing to derive and `NULL`
says so. A multivariate parent's are log-Cholesky coordinates or the
like, which nobody interprets; what the prior is about is the standard
deviations and the correlations of the effects within a block, and the
parent declares them.

## See also

[`penalty_readable()`](https://statmodels7.github.io/penalties7/reference/penalty_readable.md)
for the generic,
[`distributions7::mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.html)
for the declaration this passes through.

## Examples

``` r
# A univariate parent: the rate is the quantity, so nothing is derived.
penalty_readable(lasso_penalty(n_coef = 3), list(lambda = 1))
#> NULL

# A bivariate Gaussian prior: two standard deviations and a correlation.
mv <- distrib_penalty(
  distributions7::fixed(distributions7::mvgaussian_distrib(2),
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
