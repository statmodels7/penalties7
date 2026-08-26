# The Pieces a Marginal Criterion Consumes

The four quantities a REML or marginal likelihood criterion needs from a
quadratic penalty. `penalty_matrix()` returns the matrix \\S(\theta)\\
of the quadratic form, `penalty_rank()` its rank, `penalty_null_basis()`
an orthonormal basis of the directions it does not penalize, and
`penalty_logpdet()` the log pseudo-determinant with its first two
derivatives in the hyperparameters.

## Usage

``` r
penalty_matrix(pen, theta, ...)

penalty_rank(pen, ...)

penalty_null_basis(pen, ...)

penalty_logpdet(pen, theta, ...)
```

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object. Must be quadratic, or additive for the three it supplies.

- theta:

  A named list of hyperparameter values, or a named numeric vector
  carrying the same. `penalty_matrix()` and `penalty_logpdet()` only.

- ...:

  Passed to methods. No shipped method reads it.

## Value

`penalty_matrix()` a symmetric `q x q` matrix, where `q` is
`pen@n_coef`. `penalty_rank()` a single integer, fixed at construction.
`penalty_null_basis()` a `q x (q - r)` matrix with orthonormal columns
spanning the null space, and a `q x 0` matrix when the penalty is full
rank. The two matrices are base matrices ordinarily, and `dgCMatrix`
objects for a penalty built by
[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
with `blocks > 1`. `penalty_logpdet()` a list of three: `value`, a
single number; `grad`, a named list of one number per hyperparameter;
and `hess`, a named list of one number per hyperparameter pair, keyed
diagonals first as
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
is.

## Why the criterion needs them

A marginal criterion integrates the coefficients out under the prior the
penalty describes, and a Laplace expansion at the penalized mode leaves
\\\tfrac{1}{2}\log^{+}\lvert S(\theta)\rvert - \tfrac{1}{2}\log\lvert
H + S(\theta)\rvert\\. The first term is what these generics supply. It
is a **pseudo**-determinant, taken over the range of \\S\\ alone,
because a smoothing penalty is rank deficient by design and the ordinary
determinant is zero.

The rank and the null basis are fixed at construction and do not move
with `theta`, which is why neither takes it. That is deliberate and it
matters: counting eigenvalues of the assembled \\S(\theta)\\ above a
tolerance gives a count that falls as several smoothing parameters
spread apart, while the null space of a sum of positive semidefinite
matrices is the intersection of theirs and does not move at all.

## The log pseudo-determinant

For the plain quadratic branch \\S = \lambda D'PD\\ and

\$\$\log^{+}\lvert S \rvert = r\log\lambda + \log^{+}\lvert P
\rvert,\$\$

so the gradient is \\r/\lambda\\ and the second derivative
\\-r/\lambda^2\\, both exact and both returned. For the structured
branch the derivatives come from the parameters7 parametrization's own
log-determinant contract. Both are on the parameter scale; there is no
`scale` argument here.

## What rejects

All four reject for a penalty that is not quadratic, with a message
naming
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md).
The exception is
[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md),
which answers
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
as `FALSE` and supplies the first, second and fourth of these anyway;
`penalty_null_basis()` is the one it rejects.

## Methods

All four methods on the base class
[`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
signal an error naming
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md),
which is the predicate to ask first: a penalty that is not quadratic has
no matrix, no rank, no null basis and no log pseudo-determinant to
report, and an error is the answer to a question that does not apply.

## See also

[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
for the predicate that gates these,
[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md),
[`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md)
and
[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
for the branches that supply them,
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
for the keying of `penalty_logpdet()`'s `hess`.

## Examples

``` r
# Second differences over four coefficients: rank 2, a two-dimensional
# null space holding the level and the slope.
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
```
