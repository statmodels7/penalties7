# Construct a Sum of Quadratic Penalties

Builds the Gaussian prior whose precision is a weighted sum of fixed
matrices, one smoothing parameter to each. This is what an anisotropic
tensor-product smooth needs: one parameter per margin, so a fit may be
rough in one direction and smooth in another, where a single parameter
would force every direction to be smoothed alike.

## Usage

``` r
additive_penalty(
  mats,
  map = NULL,
  link_lambda = linkfunctions7::log_link(),
  tol = 1e-10
)
```

## Arguments

- mats:

  A list of symmetric positive semidefinite matrices, all of the same
  side, and at least one. Each is symmetrized, and a component that is
  zero, negative definite, or not symmetric to a relative `1e-8` is
  rejected by number.

- map:

  The matrix \\D\\, with as many rows as the components and one column
  per coefficient, or `NULL` (the default) for the identity. A map is
  absorbed at construction, each component becoming \\D'P_kD\\, so the
  stored `map` is always `NULL`.

- link_lambda:

  The linkfunctions7 link carrying each smoothing parameter onto the
  whole real line.
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  by default, and the same link is used for all of them.

- tol:

  The relative eigenvalue tolerance, `1e-10` by default. Used twice: to
  reject a component whose smallest eigenvalue is below `-tol` times its
  largest, and to count the rank of the normalized stack.

## Value

An
[`AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/AdditivePenalty.md)
object with one hyperparameter per component, named `lambda1`,
`lambda2`, ..., each bounded on \\(0, \infty)\\.

## The value

\$\$\rho(\beta;\lambda) = \tfrac{1}{2}\beta^\top S(\lambda)\beta -
\tfrac{1}{2}\log\mathrm{pdet}\\S(\lambda) + \tfrac{r}{2}\log 2\pi,
\qquad S(\lambda) = \sum\_{k} \lambda_k P_k,\$\$

with \\r\\ the rank of the sum. The hyperparameters are named `lambda1`,
`lambda2` and so on, in the order the components were given, and each is
positive on a log link.

## The determinant's derivatives

One eigendecomposition of \\S(\lambda)\\ gives the pseudo-inverse
\\S^{+}\\, and the quantities a marginal criterion reads follow from it:

\$\$\frac{\partial}{\partial\lambda_k}\log\mathrm{pdet}\\S =
\operatorname{tr}(S^{+}P_k), \qquad
\frac{\partial^{2}}{\partial\lambda_k\partial\lambda_l}
\log\mathrm{pdet}\\S = -\operatorname{tr}(S^{+}P_kS^{+}P_l).\$\$

Both are exact, and both agree with the assembled quantity to 0.

## The rank is not read off the sum

Counting the eigenvalues of \\S(\lambda)\\ above a tolerance gives the
right answer only while the parameters are comparable. Once they differ
by orders of magnitude the small contributions sink below any fixed
tolerance and are counted as zeros, so the rank appears to fall as the
fit is smoothed. Measured on first and second differences over five
coefficients, whose intersected null space is the constants:

|                         |                           |             |
|-------------------------|---------------------------|-------------|
| \\\lambda_2/\lambda_1\\ | eigenvalues above `1e-10` | stored rank |
| 1 to 1e8                | 4                         | 4           |
| 1e10 to 1e14            | **3**                     | 4           |

while the null basis annihilates the assembled sum to `1e-16` relative
at every one of those settings. So the rank is fixed once at
construction, from the components stacked and individually normalized,
and the object's answer cannot move. A test pins both halves.

## What this branch supplies

[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
answers `TRUE`, a sum of quadratic forms being a quadratic form, so
[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md),
[`penalty_rank()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md),
[`penalty_null_basis()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
and
[`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
all answer and
[`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
runs its three quadratic rows here. What is particular is that the
matrix moves with one hyperparameter per component, where the plain
quadratic branch has a single scale, so the log pseudo-determinant is
linear in none of them and the row that checks it compares a gradient
against `numDeriv` where the plain quadratic branch reads a slope.

[`has_prox()`](https://statmodels7.github.io/penalties7/reference/has_prox.md)
answers `FALSE`, this branch registering no
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md),
which is asked before
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
is.

## See also

[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
for one matrix with one scale,
[`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md)
for a matrix whose entries move with the hyperparameters,
[`penalty_logpdet.AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.AdditivePenalty.md)
for the determinant and its derivatives,
[`additive_sum()`](https://statmodels7.github.io/penalties7/reference/additive_sum.md)
for the decomposition they share.

## Examples

``` r
# An anisotropic tensor smooth over a 4 x 4 grid: curvature along each
# margin, penalized separately.
P <- crossprod(diff(diag(4), differences = 2))
pen <- additive_penalty(list(kronecker(diag(4), P), kronecker(P, diag(4))))
pen@params
#> [1] "lambda1" "lambda2"
c(n_coef = pen@n_coef, rank = penalty_rank(pen))
#> n_coef   rank 
#>     16     12 

# Raising one parameter charges for roughness along that margin alone.
set.seed(2)
b <- rnorm(16)
penalty_value(pen, b, list(lambda1 = 1, lambda2 = 1))
#> [1] 66.29988
penalty_value(pen, b, list(lambda1 = 1, lambda2 = 100))
#> [1] 2292.209

# With one component the branch is the plain quadratic penalty.
P2 <- crossprod(diff(diag(5), differences = 2))
bb <- c(0.4, -1.1, 0.7, 0.2, -0.3)
penalty_value(additive_penalty(list(P2)), bb, list(lambda1 = 3)) -
  penalty_value(quadratic_penalty(P2), bb, list(lambda = 3))
#> [1] 7.105427e-15

# The stored rank does not move as the parameters spread apart, where an
# eigenvalue count of the assembled sum does.
a <- additive_penalty(list(crossprod(diff(diag(5))), P2))
penalty_rank(a)
#> [1] 4
sapply(c(1, 1e10), function(r) {
  S <- penalty_matrix(a, list(lambda1 = 1, lambda2 = r))
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  sum(ev > 1e-10 * max(ev))
})
#> [1] 4 3
```
