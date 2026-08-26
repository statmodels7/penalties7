# Smoothness and Properness of an Additive Penalty

[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
returns `numeric(0)`, the value being a quadratic form and smooth
everywhere.
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
returns `TRUE` when the stored rank equals the number of coefficients.

## Arguments

- pen:

  An
  [`AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/AdditivePenalty.md)
  object.

- theta:

  A named list of `lambda1`, `lambda2`, ... Read by
  [`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md),
  whose answer does not depend on it.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
a numeric vector of length zero.
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
a single logical.

## Details

The branch registers no
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
method, so it inherits `FALSE` from
[`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
even though the value is a quadratic form. See
[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
for what that costs.

Properness is the rank of the sum against the number of coefficients,
and the rank is fixed at construction. The two components of an
anisotropic tensor smooth intersect in a null space of their own, so
such a penalty is improper: a 4 by 4 grid penalized by curvature along
each margin has rank 12 out of 16.

## See also

[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
and
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
for the generics,
[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
for why
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
answers `FALSE` here.

## Examples

``` r
P <- crossprod(diff(diag(4), differences = 2))
te <- additive_penalty(list(kronecker(diag(4), P), kronecker(P, diag(4))))

penalty_kinks(te, list(lambda1 = 1, lambda2 = 1))
#> numeric(0)
is_proper(te)
#> [1] FALSE
c(rank = penalty_rank(te), n_coef = te@n_coef)
#>   rank n_coef 
#>     12     16 

# A sum of full-rank components is proper.
is_proper(additive_penalty(list(diag(3), diag(c(2, 1, 1)))))
#> [1] TRUE
```
