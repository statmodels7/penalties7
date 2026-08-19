# Evaluate a Smoother or One of Its Derivatives

\\d^k s/du^k\\ at the given points, for \\k\\ from 0 to 5. The width is
the smoother's own unless one is supplied, which is how a consumer that
resolved the width at build evaluates without mutating the object; a
vector width is one value per point.

## Usage

``` r
smoother_deriv(smoother, u, width = NULL, order = 0L)
```

## Arguments

- smoother:

  An
  [`abs_smoother`](https://statmodels7.github.io/penalties7/reference/abs_smoother.md).

- u:

  A numeric vector.

- width:

  The width, or `NULL` for the smoother's own.

- order:

  0 to 5.

## Value

A numeric vector as long as `u`.

## See also

[`abs_smoother`](https://statmodels7.github.io/penalties7/reference/abs_smoother.md)

## Examples

``` r
smoother_deriv(smooth_hyperbolic(c = 1), 0:3, order = 1)
#> [1] 0.0000000 0.7071068 0.8944272 0.9486833
```
