# Quantities of a Quadratic Penalty for a Marginal Criterion

These functions return the four quantities that a REML or marginal
likelihood criterion requires from a quadratic penalty.
`penalty_matrix()` returns the matrix \\S(\theta)\\ of the quadratic
form, `penalty_rank()` its rank, `penalty_null_basis()` an orthonormal
basis of the directions that the penalty leaves unpenalized, and
`penalty_logpdet()` the log pseudo-determinant with its first two
derivatives in the hyperparameters.

## Usage

``` r
penalty_matrix(penalty, theta, ...)

penalty_rank(penalty, ...)

penalty_null_basis(penalty, ...)

penalty_logpdet(penalty, theta, ...)
```

## Arguments

- penalty:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object for which
  [`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
  is `TRUE`, or a separable penalty for `penalty_null_basis()`.

- theta:

  A named list of hyperparameter values, or a named numeric vector
  carrying the same. Only `penalty_matrix()` and `penalty_logpdet()`
  take this argument.

- ...:

  Passed to methods. The shipped methods do not read it.

## Value

`penalty_matrix()` returns a symmetric `q x q` matrix, where `q` is
`penalty@n_coef`. `penalty_rank()` returns a single integer, fixed at
construction. `penalty_null_basis()` returns a `q x (q - r)` matrix with
orthonormal columns spanning the null space, and a `q x 0` matrix when
there is none. The two matrices are base matrices ordinarily, and
`dgCMatrix` objects for a penalty built by
[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
with `blocks > 1`. `penalty_logpdet()` returns a list of three: `value`,
a single number; `grad`, a named list of one number per hyperparameter;
and `hess`, a named list of one number per hyperparameter pair, keyed as
in
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md),
with the diagonal pairs first.

## Role in a marginal criterion

A marginal criterion integrates the coefficients out under the prior
that the penalty defines, and a Laplace expansion at the penalized mode
leaves \\\tfrac{1}{2}\log^{+}\lvert S(\theta)\rvert -
\tfrac{1}{2}\log\lvert H + S(\theta)\rvert\\. These generics supply the
first term. It is a pseudo-determinant, taken only over the range of
\\S\\, because a smoothing penalty is rank deficient by design and its
ordinary determinant is zero.

The rank and the null basis are fixed at construction and do not move
with `theta`, so `penalty_rank()` and `penalty_null_basis()` do not take
it. A count of the eigenvalues of the assembled \\S(\theta)\\ above a
tolerance would fall as several smoothing parameters spread apart,
whereas the null space of a sum of positive semidefinite matrices is the
intersection of their null spaces and does not depend on the weights.

## The log pseudo-determinant

For the plain quadratic branch \\S = \lambda D^\top P D\\ and

\$\$\log^{+}\lvert S \rvert = r\log\lambda + \log^{+}\lvert D^\top P D
\rvert,\$\$

with \\r = \operatorname{rank}(D^\top P D)\\, so the gradient is
\\r/\lambda\\ and the second derivative \\-r/\lambda^2\\, both exact and
both returned. For the structured branch the value and the derivatives
are those of the log-determinant of the parameters7 matrix parameter
([`parameters7::param_logdet()`](https://statmodels7.github.io/parameters7/reference/param_logdet.html)
and its derivatives). All are on the parameter scale, and there is no
`scale` argument.

## Branches that supply them

The quadratic, structured and additive branches, for which
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
is `TRUE`, supply all four. For a separable penalty
([`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md)),
`penalty_null_basis()` returns an orthonormal basis of the null space of
the map, the directions of \\\beta\\ on which the penalty does not
depend, and the other three signal an error. For SCAD, MCP and any other
branch without a method of its own, all four signal an error with a
message that names
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md).

## Methods

The four methods on the base class
[`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
signal an error with a message that names
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md),
because a penalty that is not quadratic has no matrix, rank, null basis
or log pseudo-determinant. The quadratic branches register their own
methods, documented at
[`penalty_matrix.QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.QuadraticPenalty.md),
[`penalty_matrix.StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.StructuredPenalty.md)
and
[`penalty_matrix.AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.AdditivePenalty.md),
and a separable penalty registers
[`penalty_null_basis.DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_null_basis.DistribPenalty.md).

## See also

[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
for the predicate that determines which branches supply these,
[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md),
[`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md)
and
[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
for those branches,
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
for the names of `penalty_logpdet()`'s `hess`.

## Examples

``` r
# Second differences over four coefficients give rank 2 and a
# two-dimensional null space that holds the level and the slope.
pen <- quadratic_penalty(crossprod(diff(diag(4), differences = 2)))
penalty_rank(pen)
#> [1] 2
N <- penalty_null_basis(pen)
dim(N)
#> [1] 4 2

# The matrix annihilates that basis exactly.
max(abs(penalty_matrix(pen, list(lambda = 2)) %*% N))
#> [1] 1.332268e-15

# And a straight line, which the null space spans, costs nothing.
penalty_gradient(pen, c(1, 2, 3, 4), list(lambda = 2))
#> [1] 0 0 0 0

# The log pseudo-determinant is r log(lambda) plus a constant, so its
# derivatives are r/lambda and -r/lambda^2.
lp <- penalty_logpdet(pen, list(lambda = 2))
c(lp$grad$lambda, lp$hess$lambda_lambda)
#> [1]  1.0 -0.5
c(penalty_rank(pen) / 2, -penalty_rank(pen) / 4)
#> [1]  1.0 -0.5

# A separable penalty under a second-difference map has a null basis,
# and no matrix.
lsep <- lasso_penalty(map = diff(diag(4), differences = 2))
dim(penalty_null_basis(lsep))
#> [1] 4 2
try(penalty_matrix(lsep, list(lambda = 1)))
#> Error : Only a quadratic penalty exposes its matrix; see is_quadratic().
```
