# Is a Penalty a Quadratic Form?

Returns `TRUE` for
[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
(and
[`ridge_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md),
which is built on it),
[`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md)
and
[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md),
the branches whose value is the quadratic form \\\tfrac{1}{2}\beta^\top
S(\theta)\beta\\ plus a normalizing constant, with a matrix
\\S(\theta)\\ that
[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
returns. Returns `FALSE` for every other branch. The predicate
determines how a marginal likelihood criterion treats the penalty: where
it is `TRUE` the penalty's contribution to the criterion is a log
pseudo-determinant, and
[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md),
[`penalty_rank()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md),
[`penalty_null_basis()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
and
[`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
supply the pieces. Where it is `FALSE`,
[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md),
[`penalty_rank()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
and
[`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
signal an error, and
[`penalty_null_basis()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
returns a value only for a separable penalty.

## Usage

``` r
is_quadratic(penalty, ...)
```

## Arguments

- penalty:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object of any branch.

- ...:

  Passed to methods. The shipped methods do not read it.

## Value

A single logical.

## The quadratic form

A quadratic penalty carries a fixed matrix \\P\\, a map \\D\\ and one
smoothing parameter \\\lambda\\, and is the negative log-density of the
Gaussian prior on \\\beta\\ with precision \\\lambda D^\top P D\\:

\$\$\rho(\beta; \lambda) = \tfrac{\lambda}{2} (D\beta)^\top P (D\beta) -
\tfrac{1}{2}\log^{+}\lvert \lambda D^\top P D \rvert +
\tfrac{r}{2}\log(2\pi), \qquad \log^{+}\lvert \lambda D^\top P D \rvert
= r \log \lambda + \log^{+}\lvert D^\top P D \rvert,\$\$

with \\\log^{+}\\ the log pseudo-determinant and \\r =
\operatorname{rank}(D^\top P D)\\. Without a map \\D\\ is the identity.
A structured penalty has the same form with a parameters7 matrix
parameter \\S(\theta)\\ in place of \\\lambda D^\top P D\\, so its
matrix moves with several hyperparameters and its log-determinant is not
linear in any one of them.

An additive penalty has the matrix \\\sum_k \lambda_k P_k\\. A sum of
quadratic forms is a quadratic form, so it has all four quantities. Its
matrix moves with one hyperparameter per component, where the plain
quadratic branch has a single scale \\\lambda\\.

## Methods

The method on the base class
[`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
returns `FALSE`. A quadratic branch registers its own method returning
`TRUE`, as
[`is_quadratic.QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.QuadraticPenalty.md),
[`is_quadratic.StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.StructuredPenalty.md)
and
[`is_quadratic.AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.AdditivePenalty.md)
do.

## See also

[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
and the three quantities beside it,
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
and
[`has_prox()`](https://statmodels7.github.io/penalties7/reference/has_prox.md)
for the other two predicates, and
[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md),
[`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md)
and
[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
for the branches.

## Examples

``` r
# The quadratic and structured branches return TRUE.
is_quadratic(quadratic_penalty(diag(2)))
#> [1] TRUE
is_quadratic(ridge_penalty())
#> [1] TRUE
is_quadratic(structured_penalty(
  parameters7::log_cholesky(2)))
#> [1] TRUE

# The additive branch is also quadratic, and its matrix moves with one
# hyperparameter per component.
add <- additive_penalty(list(diag(3), diag(c(1, 1, 0))))
c(quadratic = is_quadratic(add), rank = penalty_rank(add))
#> quadratic      rank 
#>         1         3 

# A separable penalty, SCAD and MCP are not quadratic.
is_quadratic(lasso_penalty())
#> [1] FALSE
```
