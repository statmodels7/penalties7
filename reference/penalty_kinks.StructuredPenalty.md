# Smoothness and Kind of a Structured Penalty

[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
returns `numeric(0)`, the value being a quadratic form and smooth
everywhere.
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
returns `TRUE`, so the four marginal quantities are available.
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
returns `TRUE` when the structure has full rank.

## Arguments

- pen:

  A
  [`StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/StructuredPenalty.md)
  object.

- theta:

  A named list of the structure's free values. Read by
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
a single logical, `TRUE` when the structure's rank equals its dimension.

## Details

Properness is the structure's rank against its dimension and nothing
else. A full-rank structure is a proper Gaussian prior; a rank-deficient
one, which is admitted only as a precision, is improper and its constant
is the log pseudo-determinant over the range.

## See also

[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md),
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
and
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
for the generics,
[`penalty_matrix.StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.StructuredPenalty.md)
for the quantities
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
gates.

## Examples

``` r
pen <- structured_penalty(parameters7::ar1(4))
penalty_kinks(pen, list(log_scale = 0.2, z_rho = 0.5))
#> numeric(0)
is_quadratic(pen)
#> [1] TRUE
is_proper(pen)
#> [1] TRUE
c(penalty_rank(pen), pen@n_coef)
#> [1] 4 4
```
