# The Smallest Width a Consumer May Use

The floor of the width, derived from the expression that binds rather
than chosen. The derivatives of a smoother scale as \\s^{(k)} \sim
h^{1-k}\\, so the Jacobian column of a smoothed break-point carries
\\s''(0)/2 \sim 1/h\\ against covariate columns of order the range
\\D\\: holding the design's condition number below \\\epsilon^{-1/2}\\,
which leaves half the digits of a double to a QR of the design – the
same bound the working schedule's scaling floor rests on – gives \\h
\geq \sqrt{\epsilon}\\D\\, carried onto the width parameter's own scale
for a smoother parametrized by a squared length.

## Usage

``` r
smoother_width_floor(smoother, scale)
```

## Arguments

- smoother:

  An
  [`abs_smoother`](https://statmodels7.github.io/penalties7/reference/abs_smoother.md).

- scale:

  The scale of the covariate (its range).

## Value

A single number.

## See also

[`smoother_width`](https://statmodels7.github.io/penalties7/reference/smoother_width.md)

## Examples

``` r
smoother_width_floor(smooth_probit(), scale = 10)
#> [1] 1.490116e-07
```
