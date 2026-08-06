# Is a Penalty Quadratic?

`TRUE` only for
[`quadratic_penalty`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md),
whose matrix, rank, null basis and log pseudo-determinant the
marginal-likelihood generics expose.

## Usage

``` r
is_quadratic(pen, ...)
```

## Arguments

- pen:

  A
  [`penalty`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- ...:

  Passed to methods.

## Value

A single logical.

## Examples

``` r
is_quadratic(quadratic_penalty(diag(2)))
#> [1] TRUE
is_quadratic(ridge_penalty())
#> [1] FALSE
```
