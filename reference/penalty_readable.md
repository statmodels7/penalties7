# What a Penalty's Hyperparameters Are About

Returns the quantities a reader reads, for a penalty whose
hyperparameters are coordinates of a chart and not the quantities
themselves, together with the Jacobian from those coordinates and the
scale each one's interval belongs on. Returns `NULL` where the
hyperparameters already are the quantities, which is every branch but
one.

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

  Passed to methods. No shipped method reads it.

## Value

`NULL`, or a list of four as
[`distributions7::mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.html)
returns them: `value`, a named numeric vector of the reported
quantities; `jacobian`, their derivatives with respect to the
hyperparameters, one row per quantity and one column per hyperparameter;
`transform`, a character vector naming the scale each interval is built
on (`"log"`, `"atanh"`, `"identity"`); and `block`, a character vector
grouping the quantities for printing.

## The case this exists for

A penalty whose prior is a multivariate family carries the free values
of a matrix parameter as its hyperparameters: the logarithms of the
diagonal of a Cholesky factor and the entries below it. Nobody reads
those. What the prior is about is the standard deviations and the
correlations of the effects it describes, so
[`distributions7::mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.html)
declares them and this generic passes them through. It is the same
distinction
[`parameters7::param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.html)
makes for a matrix parameter.

The Jacobian is what a delta-method standard error needs: with \\V\\ the
variance matrix of the hyperparameters, the reported quantities have
variance \\JVJ^\top\\. The `transform` element says which scale each
quantity's confidence interval should be built on and mapped back from,
so that a standard deviation stays positive and a correlation stays
inside \\(-1, 1)\\.

## The base method

`NULL` says that the hyperparameters are the quantities and a consumer
should report them as they stand. That is the answer for every other
branch: a smoothing parameter, a rate and a shape are each read on their
own scale already.

## Methods

The method on the base class
[`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
returns `NULL`, which is the answer for every branch whose
hyperparameters already are the quantities a reader reads. Only a branch
that carries them on a chart registers a method.

## See also

[`penalty_value()`](https://statmodels7.github.io/penalties7/reference/penalty_value.md),
[`distributions7::mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.html)
for the declaration this passes through,
[`parameters7::param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.html)
for the same distinction on a matrix parameter.

## Examples

``` r
# A bivariate Gaussian prior on three pairs of effects. Its hyperparameters
# are log-Cholesky coordinates, which nobody interprets.
pen <- distrib_penalty(
  distributions7::fixed(distributions7::mvgaussian1_distrib(2),
                        mu1 = 0, mu2 = 0), n_coef = 6)
pen@params
#> [1] "sigma_log_L1" "sigma_log_L2" "sigma_L2.1"  

# What the prior is about is two standard deviations and a correlation.
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

# A smoothing parameter is already the quantity it names, so there is
# nothing to derive.
penalty_readable(quadratic_penalty(diag(2)), list(lambda = 1))
#> NULL
```
