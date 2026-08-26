# Print a Penalty

Writes one line naming the penalty, the number of coefficients it takes,
the number of rows its map returns, and its hyperparameters. The row
count is `nrow(map)`, or `n_coef` when the map is `NULL`, so an unmapped
penalty shows the same number twice.

## Arguments

- x:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object of any branch.

- ...:

  Unused, accepted for consistency with
  [`print()`](https://rdrr.io/r/base/print.html).

## Value

`x`, invisibly.

## Examples

``` r
quadratic_penalty(diag(3))
#> quadratic penalty on 3 coefficient(s) through 3 row(s); theta: lambda

# A map narrows what the penalty sees: three coefficients, two differences.
quadratic_penalty(diag(2), map = diff(diag(3)))
#> quadratic penalty on 3 coefficient(s) through 2 row(s); theta: lambda

# A penalty with no hyperparameters says so.
distrib_penalty(
  distributions7::fixed(distributions7::gaussian1_distrib(),
                        mu = 0, sigma = 1),
  n_coef = 3)
#> separable [fixed gaussian1 [mu=0,sigma=1]] penalty on 3 coefficient(s) through 3 row(s); theta: (none)
```
