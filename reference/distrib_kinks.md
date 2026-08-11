# Where the Parent's Log-Density Has a Kink

The points a separable penalty is not differentiable at, derived from
the distribution it is built on.

## Usage

``` r
distrib_kinks(d)
```

## Arguments

- d:

  A univariate distributions7 object.

## Value

A numeric vector, possibly empty.

## Details

A distribution records which of its parameters the log-likelihood is
differentiable in, through `params_smooth`, and for a location family a
location that is not smooth is a kink in the argument at that location.
A penalty is the negative log-density read in the coefficient, so a
[`fixed`](https://statmodels7.github.io/distributions7/reference/fixed.html)
wrapper holding such a parameter at a value puts the kink there:
`fixed(laplace_distrib(), mu = 0)` is the lasso and has a kink at zero.

Each candidate is the value its parameter is held at, and whether it is
a kink is then measured rather than inferred, by comparing the one-sided
derivatives of the log-density across it. Inferring alone would put a
kink on any family whose non-smooth parameter is not a location;
measuring alone would need somewhere to look. A candidate whose
derivative does not jump is dropped.

A parent that declares every parameter smooth, or that fixes none of the
ones it declares non-smooth, has no candidate and gets `numeric(0)`.
Nothing is taken from a parameter that is free, its value being whatever
the hyperparameters say at the time.

## See also

[`distrib_penalty`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md),
[`penalty_kinks`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)

## Examples

``` r
distrib_kinks(distributions7::fixed(distributions7::laplace_distrib(),
                                    mu = 0))
#> [1] 0
distrib_kinks(distributions7::fixed(distributions7::gaussian1_distrib(),
                                    mu = 0))
#> numeric(0)
```
