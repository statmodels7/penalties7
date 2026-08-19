# Check a Smoother of the Absolute Value Numerically

The sibling of `check_penalty` for user-written smoothers: the
structural properties of \\s\\ (even, with an odd first derivative
bounded by one, convex, and matching \\\lvert u\rvert\\ in the tails)
and each derivative order against one numerical differentiation of the
analytic order below it, never a nested difference.

## Usage

``` r
check_abs_smoother(smoother, width = NULL, tol = 1e-06, verbose = TRUE)
```

## Arguments

- smoother:

  An
  [`abs_smoother`](https://statmodels7.github.io/penalties7/reference/abs_smoother.md).

- width:

  The width to check at, or `NULL` for the smoother's own (0.5 where it
  carries none).

- tol:

  The comparison tolerance.

- verbose:

  Logical; print the table.

## Value

A data frame with one row per check, invisibly when printed.

## Details

The grid covers the transition and the tails, placed by the smoother's
own intrinsic scale \\s(0)\\. Order \\k\\ is compared against numDeriv
applied ONCE to the analytic order \\k - 1\\, so a wrong derivative is
caught against an independent route while the reference never
degenerates into a difference of differences.

## See also

[`abs_smoother`](https://statmodels7.github.io/penalties7/reference/abs_smoother.md),
[`check_penalty`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)

## Examples

``` r
res <- check_abs_smoother(smooth_probit(h = 0.3))
#> check_abs_smoother: probit at h = 0.3
#>                                         check    max_error status
#>                                     s is even 6.023629e-17     OK
#>                                     s' is odd 1.110223e-16     OK
#>                           |s'| bounded by one 0.000000e+00     OK
#>                                      s convex 0.000000e+00     OK
#>                      matches |u| in the tails 0.000000e+00     OK
#>                order 1 vs numDeriv on order 0 5.906771e-11     OK
#>                order 2 vs numDeriv on order 1 6.359873e-11     OK
#>                order 3 vs numDeriv on order 2 1.267688e-10     OK
#>                order 4 vs numDeriv on order 3 1.657041e-11     OK
#>                order 5 vs numDeriv on order 4 4.669206e-11     OK
#>                 tau_correction bounded by tau 0.000000e+00     OK
#>  tau_correction is the identity at width zero 0.000000e+00     OK
all(res$status == "OK")
#> [1] TRUE
```
