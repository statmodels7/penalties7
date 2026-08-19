# The Probit Smoother of the Absolute Value

\\s(u) = \mathbb{E}\lvert u + hZ\rvert\\ for \\Z\\ standard normal:

\$\$s(u) = u\\\bigl(2\Phi(u/h) - 1\bigr) + 2h\\\phi(u/h), \qquad s'(u) =
2\Phi(u/h) - 1, \qquad s''(u) = \frac{2}{h}\\\phi(u/h),\$\$

the recommended default. Its excess over \\\lvert u\rvert\\ has gaussian
tails, so the smoothing bias is confined to a window of width \\h\\
around the kink, and it satisfies an exact convolution identity:
smoothing the step with width \\h\\ is convolving the break-point with
\\N(0, h^2)\\, so the apparent scale of a random break-point composes as
\\\tau^2\_{\mathrm{apparent}} = \tau^2 + h^2\\ and the smoother declares
the closed correction \\\tau\_{\mathrm{true}} = \sqrt{\tau^2 - h^2}\\.

## Usage

``` r
smooth_probit(h = NULL, per_group = FALSE)
```

## Arguments

- h:

  The transition width, or `NULL` (the default) to be resolved at build
  from the covariate's spacing.

- per_group:

  Whether the width is resolved per group.

## Value

An
[`abs_smoother`](https://statmodels7.github.io/penalties7/reference/abs_smoother.md).

## See also

[`abs_smoother`](https://statmodels7.github.io/penalties7/reference/abs_smoother.md),
[`smooth_hyperbolic`](https://statmodels7.github.io/penalties7/reference/smooth_hyperbolic.md),
[`smooth_quintic`](https://statmodels7.github.io/penalties7/reference/smooth_quintic.md)

## Examples

``` r
sm <- smooth_probit(h = 0.2)
smoother_deriv(sm, 0, order = 0)  # 2 * h * dnorm(0)
#> [1] 0.1595769
```
