# What a Penalty's Hyperparameters Are About

The quantities a reader reads, where the hyperparameters are coordinates
of a chart rather than the quantities themselves, with the Jacobian from
those coordinates and the scale each one's interval belongs on.

## Usage

``` r
penalty_readable(pen, theta, ...)
```

## Arguments

- pen:

  A
  [`penalty`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- theta:

  A named list of hyperparameter values.

- ...:

  Passed to methods.

## Value

`NULL`, or a list with `value`, `jacobian`, `transform` and `block`, as
[`mv_derived`](https://statmodels7.github.io/distributions7/reference/mv_derived.html)
returns them.

## Details

The case this exists for is a penalty whose prior is a multivariate
family: its hyperparameters are the free values of a matrix parameter –
the logarithms of the diagonal of a Cholesky factor and the entries
below it – and nobody reads those. What the prior is about is the
standard deviations and the correlations of the effects it describes,
and
[`mv_derived`](https://statmodels7.github.io/distributions7/reference/mv_derived.html)
declares them, so this is the same distinction
[`param_readable`](https://statmodels7.github.io/parameters7/reference/param_readable.html)
makes for a matrix parameter and `term_readable` for a fitted term.

The base method returns `NULL`, which says that the hyperparameters ARE
the quantities and a consumer should report them as they stand. That is
the honest answer for every other branch: a smoothing parameter, a rate,
a shape are each read on their own scale already.

## See also

[`penalty_value`](https://statmodels7.github.io/penalties7/reference/penalty_value.md),
[`mv_derived`](https://statmodels7.github.io/distributions7/reference/mv_derived.html)

## Examples

``` r
pen <- distrib_penalty(
  distributions7::fixed(distributions7::mvgaussian_distrib(2),
                        mu1 = 0, mu2 = 0), n_coef = 6)
penalty_readable(pen, list(sigma_log_L1 = 0.2, sigma_log_L2 = -0.1,
                           sigma_L2.1 = 0.5))$value
#>     sd_v1     sd_v2 cor_v1_v2 
#> 1.2214028 1.0337943 0.4836552 

# a smoothing parameter is already the quantity it names
penalty_readable(quadratic_penalty(diag(2)), list(lambda = 1))
#> NULL
```
