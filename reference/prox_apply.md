# Apply a Piecewise Linear Table

Evaluates in R the map that
[`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
describes. A compiled coordinate descent applies the same map to each
coordinate, and the tests compare the result with
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md).

## Usage

``` r
prox_apply(spec, u)
```

## Arguments

- spec:

  A table, as
  [`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
  returns it. Its matrices must have at least as many rows as `u` is
  long.

- u:

  The points, one per coefficient, a numeric vector.

## Value

A numeric vector as long as `u`.

## Details

For coordinate \\j\\ the function finds the first cut that \\\lvert
u_j\rvert\\ does not exceed and returns
\\\mathrm{sign}(u_j)(a\_{jk}\lvert u_j\rvert + b\_{jk})\\ for that
piece. The last cut of every shipped table is `Inf`, so a piece is
always found. If \\\lvert u_j\rvert\\ exceeds every cut of its row,
which can happen only in a table whose last cut is finite, the last
piece is used.

A table is verified by comparing this function with
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
at every cut, on both sides of every cut and over a grid of points. For
SCAD and MCP at a step close to its limit the table is less accurate
than
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md),
as the page of
[`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
explains.

## See also

[`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
for the table,
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
for the operator in whose place the table is used.

## Examples

``` r
pen <- lasso_penalty(n_coef = 2L)
sp <- penalty_prox_spec(pen, list(lambda = 1.5), step = c(0.5, 2))

# The first coordinate is thresholded at 0.75 and the second at 3.
prox_apply(sp, c(2, 2))
#> [1] 1.25 0.00
prox_apply(sp, c(0.5, 2.9))
#> [1] 0 0

# This is exactly what the operator gives, coordinate by coordinate.
one <- lasso_penalty(n_coef = 1)
c(penalty_prox(one, 2, 0.5, list(lambda = 1.5)),
  penalty_prox(one, 2, 2, list(lambda = 1.5)))
#> [1] 1.25 0.00

# The map is odd, so negating the point negates the answer.
prox_apply(sp, c(-2, -2))
#> [1] -1.25  0.00
```
