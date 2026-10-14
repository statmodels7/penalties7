# Smoothness and Properness of a Separable Penalty

[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
returns the kinks recorded on the object at construction, and
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
returns `TRUE` when there is no map or the map is square and
nonsingular, and `FALSE` otherwise.

## Arguments

- pen:

  A
  [`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
  object.

- theta:

  A named list of the parent's free parameters. It is accepted by
  [`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md),
  whose result does not depend on it.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
a numeric vector, `pen@kinks` as it stands.
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
a single logical.

## Details

The kinks do not move with `theta`, unlike those of SCAD and MCP: they
are fixed by the value at which the location of the parent is held, and
they are set when the penalty is built, by
[`distrib_kinks()`](https://statmodels7.github.io/penalties7/reference/distrib_kinks.md)
or by the argument `kinks` of
[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md).
Of the named penalties, the lasso and the elastic net report `0` and
[`heavy_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md)
reports `numeric(0)`.

Properness is a test on the map. The parent is a density, so
\\\exp(-\rho)\\ integrates to one over the values at which the penalty
is read, and under a square nonsingular map the value carries
\\-\log\lvert\det D\rvert\\, so that \\\exp(-\rho)\\ integrates to one
over \\\beta\\ as well. A map whose rank is below its number of columns
leaves directions of \\\beta\\ on which the penalty does not depend, and
along those directions the integral diverges. Under a map with more rows
than columns and full column rank the integral does not equal one, and
it depends on the hyperparameters through a constant that the value does
not carry. The test reads the `log_jacobian` stored at construction, or
the entries of a diagonal map, through
[dp_log_jacobian()](https://statmodels7.github.io/penalties7/reference/dp_diag_entries.md),
and returns `TRUE` exactly when that quantity is finite.

## See also

[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
and
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
for the generics,
[`distrib_kinks()`](https://statmodels7.github.io/penalties7/reference/distrib_kinks.md)
for the derivation of the kinks,
[`penalty_null_basis.DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_null_basis.DistribPenalty.md)
for the directions that an improper map leaves free,
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
for the proximal operator at a kink.

## Examples

``` r
# The kink of the lasso is at zero for every lambda.
penalty_kinks(lasso_penalty(n_coef = 3), list(lambda = 1))
#> [1] 0
penalty_kinks(lasso_penalty(n_coef = 3), list(lambda = 20))
#> [1] 0

# A Student t prior has none.
penalty_kinks(heavy_penalty(n_coef = 3), list(sigma = 1, nu = 4))
#> numeric(0)

# Proper without a map and with a square nonsingular map; improper with a
# map that has more rows than columns, with a singular one, and with one
# that has fewer rows than columns.
c(none     = is_proper(lasso_penalty(n_coef = 3)),
  square   = is_proper(lasso_penalty(map = matrix(c(2, 1, 0, 3), 2))),
  tall     = is_proper(lasso_penalty(map = rbind(diag(2), c(1, 1)))),
  singular = is_proper(lasso_penalty(map = matrix(1, 2, 2))),
  short    = is_proper(lasso_penalty(map = diff(diag(3), differences = 2))))
#>     none   square     tall singular    short 
#>     TRUE     TRUE    FALSE    FALSE    FALSE 
```
