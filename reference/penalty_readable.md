# Interpretable Quantities of a Penalty's Hyperparameters

For a penalty whose hyperparameters are coordinates of a chart rather
than interpretable quantities, returns the quantities to be reported,
together with the Jacobian from those coordinates and the scale on which
each quantity's interval is built. Returns `NULL` where the
hyperparameters are reported as they stand.

## Usage

``` r
penalty_readable(pen, theta, ...)
```

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object of any branch.

- theta:

  A named list of hyperparameter values, or a named numeric vector
  carrying the same, holding every name in `pen@params`.

- ...:

  Passed to methods. The shipped methods do not read it.

## Value

`NULL`, or a list of four as
[`distributions7::mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.html)
returns them: `value`, a named numeric vector of the reported
quantities; `jacobian`, their derivatives with respect to the
hyperparameters, one row per quantity and one column per hyperparameter;
`transform`, a character vector naming the scale each interval is built
on (`"log"`, `"atanh"`, `"identity"`); and `block`, a character vector
grouping the quantities for printing.

## Multivariate priors

A separable penalty whose parent is a multivariate family carries the
free values of a matrix parameter as its hyperparameters: the logarithms
of the diagonal of a Cholesky factor and the entries below it, which are
not interpretable directly. The interpretable quantities are the
standard deviations and the correlations of the effects, which
[`distributions7::mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.html)
computes and this generic returns. The same distinction is made by
[`parameters7::param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.html)
for a matrix parameter.

The Jacobian gives a delta-method standard error: with \\V\\ the
variance matrix of the hyperparameters, the reported quantities have
variance \\JVJ^\top\\. The `transform` element names the scale on which
each quantity's confidence interval is built and from which it is mapped
back, so that a standard deviation stays positive and a correlation
stays inside \\(-1, 1)\\.

## The base method

`NULL` means that the hyperparameters are reported as they stand. This
is the result for the quadratic, additive, SCAD and MCP branches and for
a separable penalty with a univariate parent, whose hyperparameters (a
smoothing parameter, a rate, a scale, a shape) are interpretable on
their own scale. A structured penalty returns `NULL` when its structure
declares no quantities through
[`parameters7::param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.html),
as
[`parameters7::log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.html)
and
[`parameters7::inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.html)
do.

## Methods

The method on the base class
[`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
returns `NULL`. A separable penalty registers
[`penalty_readable.DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_readable.DistribPenalty.md),
which returns the result of
[`distributions7::mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.html)
for a multivariate parent and `NULL` for a univariate one. A structured
penalty registers
[`penalty_readable.StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_readable.StructuredPenalty.md),
which returns the quantities that
[`parameters7::param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.html)
declares for its structure, in the same shape and under names prefixed
with `precision_`, and `NULL` for a structure that declares none.

## See also

[`penalty_value()`](https://statmodels7.github.io/penalties7/reference/penalty_value.md),
[`distributions7::mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.html)
for the quantities returned for a multivariate prior,
[`parameters7::param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.html)
for the same distinction on a matrix parameter.

## Examples

``` r
# A bivariate Gaussian prior on three pairs of effects. Its hyperparameters
# are log-Cholesky coordinates.
pen <- distrib_penalty(
  distributions7::fixed(distributions7::mvgaussian1_distrib(2),
                        mu1 = 0, mu2 = 0), n_coef = 6)
pen@params
#> [1] "sigma_log_L1" "sigma_log_L2" "sigma_L2.1"  

# The reported quantities are two standard deviations and a correlation.
r <- penalty_readable(pen, list(sigma_log_L1 = 0.2, sigma_log_L2 = -0.1,
                                sigma_L2.1 = 0.5))
r$value
#>     sd_v1     sd_v2 cor_v1_v2 
#> 1.2214028 1.0337943 0.4836552 
r$transform
#>     sd_v1     sd_v2 cor_v1_v2 
#>     "log"     "log"   "atanh" 
dim(r$jacobian)
#> [1] 3 3

# A smoothing parameter is reported as it stands.
penalty_readable(quadratic_penalty(diag(2)), list(lambda = 1))
#> NULL
```
