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

## See also

[`is_quadratic`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md),
[`has_prox`](https://statmodels7.github.io/penalties7/reference/has_prox.md),
[`penalty_matrix`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md),
[`penalty_rank`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md),
[`penalty_null_basis`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md),
[`penalty_logpdet`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)

## Examples

``` r
is_proper(quadratic_penalty(diag(2)))
#> [1] TRUE
is_proper(scad_penalty())
#> [1] FALSE
```
