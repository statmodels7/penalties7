# Smoothness and Properness of a Separable Penalty

[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
returns the kinks recorded on the object at construction, and
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
returns `TRUE` unless the map has fewer rows than columns.

## Arguments

- pen:

  A
  [`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
  object.

- theta:

  A named list of the parent's free parameters. Read by
  [`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md),
  whose answer does not depend on it.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
a numeric vector, `pen@kinks` as it stands.
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
a single logical.

## Details

The kinks do not move with `theta`, unlike SCAD's and MCP's: they are
properties of where the parent's location is held, settled by
[`distrib_kinks()`](https://statmodels7.github.io/penalties7/reference/distrib_kinks.md)
when the penalty was built. The lasso and the elastic net report `0`;
every other shipped parent reports `numeric(0)`.

Properness here is a test on the map's shape. The parent is a density by
construction, so \\\exp(-\rho)\\ integrates over the values the penalty
reads; what can fail is that a map with fewer rows than columns leaves
directions of \\\beta\\ the penalty never sees, and along those the
integral diverges. With no map, or with a map at least as tall as it is
wide, the answer is `TRUE`.

## See also

[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
and
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
for the generics,
[`distrib_kinks()`](https://statmodels7.github.io/penalties7/reference/distrib_kinks.md)
for where the kinks came from,
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
for the operator a kink calls for.

## Examples

``` r
# The lasso has a kink at zero and does not move with lambda.
penalty_kinks(lasso_penalty(n_coef = 3), list(lambda = 1))
#> [1] 0
penalty_kinks(lasso_penalty(n_coef = 3), list(lambda = 20))
#> [1] 0

# A Gaussian or Student t prior has none.
penalty_kinks(heavy_penalty(n_coef = 3), list(sigma = 1, nu = 4))
#> numeric(0)

# A map that is shorter than it is wide leaves directions unpenalized.
c(short = is_proper(lasso_penalty(map = diff(diag(3), differences = 2))),
  none  = is_proper(lasso_penalty(n_coef = 3)))
#> short  none 
#> FALSE  TRUE 
```
