# The Hyperbolic Smoother of the Absolute Value

\\s(u) = \sqrt{u^2 + c}\\, the simplest smoother, whose parameter \\c\\
is a squared length (the transition width in covariate units is
\\\sqrt{c}\\). Its excess over \\\lvert u\rvert\\ decays only as
\\c/(4\lvert u\rvert)\\ in the tails, so the smoothing bias spreads away
from the kink – the worst bias profile of the three at equal width – and
no convolution identity relates \\c\\ to the scale of a random
break-point, so `tau_correction` is `NULL` and a consumer reports the
apparent scale alone.

## Usage

``` r
smooth_hyperbolic(c = NULL)
```

## Arguments

- c:

  The squared transition width, or `NULL` (the default) to be resolved
  at build as the square of the covariate's spacing.

## Value

An
[`abs_smoother`](https://statmodels7.github.io/penalties7/reference/abs_smoother.md).

## See also

[`abs_smoother`](https://statmodels7.github.io/penalties7/reference/abs_smoother.md),
[`smooth_probit`](https://statmodels7.github.io/penalties7/reference/smooth_probit.md)

## Examples

``` r
sm <- smooth_hyperbolic(c = 0.04)
smoother_deriv(sm, 0, order = 0)  # sqrt(c)
#> [1] 0.2
```
