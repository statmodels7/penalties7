# Does a Penalty Supply a Proximal Operator?

`TRUE` when
[`penalty_prox`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
can be evaluated for this penalty at a suitable step, so that a caller
may choose a proximal method over a smooth one without provoking an
error.

## Usage

``` r
has_prox(pen)
```

## Arguments

- pen:

  A
  [`penalty`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

## Value

A single logical.

## See also

[`penalty_prox`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)

## Examples

``` r
c(lasso = has_prox(lasso_penalty()), scad = has_prox(scad_penalty()))
#> lasso  scad 
#>  TRUE  TRUE 
```
