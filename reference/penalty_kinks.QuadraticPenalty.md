# Smoothness and Kind of a Quadratic Penalty

The three questions a consumer asks before choosing how to fit a
quadratic block.
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
returns `numeric(0)`, the value being a polynomial and smooth
everywhere.
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
returns `TRUE`, so the four marginal quantities are available.
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
returns `TRUE` only when the penalty has full rank.

## Arguments

- pen:

  A
  [`QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/QuadraticPenalty.md)
  object.

- theta:

  A list holding `lambda`. Read by
  [`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md),
  whose answer does not depend on it.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
a numeric vector of length zero.
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
the single logical `TRUE`.
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
a single logical, `TRUE` when the rank equals the number of coefficients
and the null basis has no columns.

## Details

Properness is a rank test and nothing else: \\\exp(-\rho)\\ is flat
along every null direction of \\D'PD\\, so it integrates only when there
are none. A ridge over \\q\\ coefficients is proper; second differences
over \\q\\ leave the level and the slope free, have rank \\q - 2\\, and
are not. Both are usable, and the difference shows in the normalizing
constant, which is taken over the range alone.

Because there are no kinks, a quadratic block goes to whatever smooth
method the fit uses, and
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
is a linear solve rather than a threshold.

## See also

[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md),
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
and
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
for the generics,
[`penalty_matrix.QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.QuadraticPenalty.md)
for the quantities
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
gates.

## Examples

``` r
ridge <- quadratic_penalty(diag(4))
curve <- quadratic_penalty(crossprod(diff(diag(4), differences = 2)))

penalty_kinks(ridge, list(lambda = 1))
#> numeric(0)
c(is_quadratic(ridge), is_quadratic(curve))
#> [1] TRUE TRUE
c(is_proper(ridge), is_proper(curve))
#> [1]  TRUE FALSE
c(penalty_rank(ridge), penalty_rank(curve))
#> [1] 4 2
```
