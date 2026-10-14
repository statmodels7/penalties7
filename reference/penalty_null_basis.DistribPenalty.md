# Null Basis of a Separable Penalty

Returns an orthonormal basis of the null space of the map \\D\\, that
is, of the directions of \\\beta\\ on which the penalty does not depend.
The basis has no columns when there is no map or when the map has full
column rank. Under a diagonal map the basis is made of the unit vectors
at the entries of the diagonal that are zero or missing, so it is
spanned by coordinates.

## Arguments

- pen:

  A
  [`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
  object.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric matrix with `pen@n_coef` rows and one column per null
direction.

## Details

A separable penalty is not a quadratic form, and its log-density is not
computed through a pseudo-determinant. A marginal criterion needs from
an improper separable prior only the directions that the prior leaves
free, and this basis supplies them. The rank of a non-diagonal map is
the rank that [`qr()`](https://rdrr.io/r/base/qr.html) reports at its
default tolerance, and the basis is formed by the right singular vectors
of the map beyond that rank. For a square map this is the same rank test
that decides
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md).
A map with more rows than columns and full column rank has an empty
basis, although
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
returns `FALSE` for it.

## See also

[`penalty_null_basis()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
for the generic,
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
for the test of properness.

## Examples

``` r
# No map, and a diagonal map with a zero at the second coordinate.
penalty_null_basis(lasso_penalty(n_coef = 3))
#>     
#> [1,]
#> [2,]
#> [3,]
penalty_null_basis(lasso_penalty(map = Matrix::Diagonal(x = c(1, 0, 2))))
#>      [,1]
#> [1,]    0
#> [2,]    1
#> [3,]    0

# A map with one row and three columns leaves two directions free.
penalty_null_basis(heavy_penalty(map = diff(diag(3), differences = 2)))
#>           [,1]       [,2]
#> [1,] 0.8164966 -0.4082483
#> [2,] 0.5265986  0.2367007
#> [3,] 0.2367007  0.8816497
```
