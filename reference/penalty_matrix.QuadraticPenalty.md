# Marginal Quantities of a Quadratic Penalty

The four pieces a REML or marginal likelihood criterion reads.
[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
returns \\\lambda D'PD\\,
[`penalty_rank()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
the rank fixed at construction,
[`penalty_null_basis()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
the stored orthonormal basis of the null space, and
[`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
the log pseudo-determinant with its first two derivatives in
\\\lambda\\.

## Arguments

- pen:

  A
  [`QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/QuadraticPenalty.md)
  object.

- theta:

  A list holding `lambda`.
  [`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
  and
  [`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
  only.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
a symmetric matrix of side `pen@n_coef`, a base matrix ordinarily and a
`dgCMatrix` when the penalty was built with `blocks > 1`.
[`penalty_rank()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
a single integer.
[`penalty_null_basis()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
a matrix with `pen@n_coef` rows and `pen@n_coef - penalty_rank(pen)`
orthonormal columns, sparse under `blocks > 1`.
[`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
a list of `value` (a number), `grad` (a list of one number, `lambda`)
and `hess` (a list of one number, `lambda_lambda`).

## Details

Since \\S(\lambda) = \lambda D'PD\\ and only the range contributes,

\$\$\log^{+}\lvert S(\lambda) \rvert = r\log\lambda +
\log\mathrm{pdet}(P), \qquad \frac{\partial}{\partial\lambda} =
\frac{r}{\lambda}, \qquad \frac{\partial^2}{\partial\lambda^2} =
-\frac{r}{\lambda^2},\$\$

all three exact and none of them needing a decomposition at call time.
The rank and the null basis do not depend on \\\lambda\\ and their
generics take no `theta`.

## See also

[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
and its three siblings for the generics,
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
for the predicate that gates them,
[`penalty_value.QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.QuadraticPenalty.md)
for the value that carries the same constant.

## Examples

``` r
pen <- quadratic_penalty(crossprod(diff(diag(5), differences = 2)))
penalty_rank(pen)
#> [1] 3

# The null basis spans the straight lines, and the matrix kills it.
N <- penalty_null_basis(pen)
dim(N)
#> [1] 5 2
max(abs(penalty_matrix(pen, list(lambda = 3)) %*% N))
#> [1] 8.382184e-15

# The log pseudo-determinant is r log(lambda) plus a constant.
lp <- penalty_logpdet(pen, list(lambda = 3))
c(lp$grad$lambda, lp$hess$lambda_lambda)
#> [1]  1.0000000 -0.3333333
c(penalty_rank(pen) / 3, -penalty_rank(pen) / 9)
#> [1]  1.0000000 -0.3333333
```
