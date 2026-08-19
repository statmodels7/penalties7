# Resolve a Smoother's Width from a Spacing

The width the smoother carries where it carries one, and otherwise the
given spacing carried onto the width parameter's own scale (the square,
for a smoother parametrized by a squared length). It is what a consumer
calls at build: a break-point term hands it the median spacing of its
covariate, the smallest transition the data can tell from a step.

## Usage

``` r
smoother_width(smoother, spacing)
```

## Arguments

- smoother:

  An
  [`abs_smoother`](https://statmodels7.github.io/penalties7/reference/abs_smoother.md).

- spacing:

  A spacing in covariate units, one value or one per group.

## Value

The width, the same length as `spacing` when resolved from it.

## See also

[`smoother_width_floor`](https://statmodels7.github.io/penalties7/reference/smoother_width_floor.md)

## Examples

``` r
smoother_width(smooth_probit(), 0.3)          # 0.3
#> [1] 0.3
smoother_width(smooth_hyperbolic(), 0.3)      # 0.09
#> [1] 0.09
smoother_width(smooth_probit(h = 0.5), 0.3)   # 0.5, the width it holds
#> [1] 0.5
```
