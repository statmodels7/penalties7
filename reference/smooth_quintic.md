# The Quintic Smoother of the Absolute Value

\\\lvert u\rvert\\ replaced inside \\\[-h, h\]\\ by the even polynomial
matching it to third order at the seam, and EXACT outside: with \\t =
u/h\\,

\$\$s(u) = \frac{h}{16}\\\bigl(5 + 15t^2 - 5t^4 + t^6\bigr) \quad
\text{for } \lvert t\rvert \< 1, \qquad s(u) = \lvert u\rvert \quad
\text{otherwise},\$\$

so \\s''(u) = 15(1 - t^2)^2/(8h)\\ vanishes to second order at \\\pm h\\
and the function is \\C^3\\, with a jump in the fourth derivative at the
seam. The smoothing bias is exactly zero outside \\\[\psi - h, \psi +
h\]\\, which is the cleanest fixed-width choice; the branches are why
the contract carries the derivatives as functions,
[`deriv`](https://rdrr.io/r/stats/deriv.html) not reading a clamp.

## Usage

``` r
smooth_quintic(h = NULL)
```

## Arguments

- h:

  The transition half-width, or `NULL` (the default) to be resolved at
  build from the covariate's spacing.

## Value

An
[`abs_smoother`](https://statmodels7.github.io/penalties7/reference/abs_smoother.md).

## See also

[`abs_smoother`](https://statmodels7.github.io/penalties7/reference/abs_smoother.md),
[`smooth_probit`](https://statmodels7.github.io/penalties7/reference/smooth_probit.md)

## Examples

``` r
sm <- smooth_quintic(h = 0.5)
smoother_deriv(sm, 1, order = 0)  # exactly 1: outside the transition
#> [1] 1
```
