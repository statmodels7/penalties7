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

## Details

The operator in question is

\$\$\operatorname{prox}\_{t\rho}(v) = \arg\min\_{\beta} \Bigl\\
\tfrac{1}{2}\lVert \beta - v \rVert^{2} + t\\ \rho(\beta; \theta)
\Bigr\\,\$\$

which a proximal gradient method evaluates once per iteration and which
therefore has to be available in closed form, or as a solve, for the
method to be worth using. A penalty carries one when it is quadratic or
structured (one linear solve at any map), when it is separable with a
parent whose operator is closed or whose stationarity condition has a
coordinatewise root, or when it is SCAD or MCP over their convex
regions.

## See also

[`penalty_prox`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)

## Examples

``` r
c(lasso = has_prox(lasso_penalty()), scad = has_prox(scad_penalty()))
#> lasso  scad 
#>  TRUE  TRUE 
```
