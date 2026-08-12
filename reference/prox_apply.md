# Apply a Piecewise Linear Table

Evaluates the map
[`penalty_prox_spec`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
describes, in R, which is what a compiled loop does and what the tests
compare the table against.

## Usage

``` r
prox_apply(spec, u)
```

## Arguments

- spec:

  A table, as
  [`penalty_prox_spec`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
  returns it.

- u:

  The points, one per coefficient.

## Value

A numeric vector as long as `u`.

## See also

[`penalty_prox_spec`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)

## Examples

``` r
pen <- lasso_penalty(n_coef = 2L)
sp <- penalty_prox_spec(pen, list(lambda = 1.5), step = c(0.5, 2))
prox_apply(sp, c(2, 2))
#> [1] 1.25 0.00
```
