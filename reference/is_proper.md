# Is a Penalty a Proper Prior?

`TRUE` when \\\exp(-\rho)\\ integrates to one over the penalized
coordinates, so that the value is exactly a negative log-density;
`FALSE` for the improper ones (a rank-deficient quadratic, SCAD, MCP),
whose value is the bare \\\rho\\.

## Usage

``` r
is_proper(pen, ...)
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
is_proper(quadratic_penalty(diag(2)))
#> [1] TRUE
is_proper(scad_penalty())
#> [1] FALSE
```
