# Apply a Piecewise Linear Table

Evaluates the map
[`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
describes, in R. This is what a compiled coordinate descent does to each
coordinate, and what the tests compare the table against.

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

For coordinate \\j\\ it finds the first cut that \\\lvert u_j\rvert\\
does not exceed and returns \\\mathrm{sign}(u_j)(a\_{jk}\lvert
u_j\rvert + b\_{jk})\\ for that piece. The last cut is `Inf`, so a piece
is always found; the `NA` guard covers a table whose cuts are not
increasing, which no shipped branch produces.

Comparing this against
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
is how a table is verified. Over the five families that have one, probed
at every cut, on both sides of every cut and across the whole range, the
worst disagreement is `8.9e-16`.

## See also

[`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
for the table,
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
for the operator it stands in for.

## Examples

``` r
pen <- lasso_penalty(n_coef = 2L)
sp <- penalty_prox_spec(pen, list(lambda = 1.5), step = c(0.5, 2))

# The first coordinate is thresholded at 0.75 and the second at 3.
prox_apply(sp, c(2, 2))
#> [1] 1.25 0.00
prox_apply(sp, c(0.5, 2.9))
#> [1] 0 0

# Which is exactly what the operator gives, coordinate by coordinate.
one <- lasso_penalty(n_coef = 1)
c(penalty_prox(one, 2, 0.5, list(lambda = 1.5)),
  penalty_prox(one, 2, 2, list(lambda = 1.5)))
#> [1] 1.25 0.00

# The map is odd, so negating the point negates the answer.
prox_apply(sp, c(-2, -2))
#> [1] -1.25  0.00
```
