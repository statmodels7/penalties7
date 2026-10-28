# Construct a Sum of Quadratic Penalties

Builds the Gaussian prior whose precision is a weighted sum of fixed
positive semidefinite matrices, with one smoothing parameter for each
matrix. An anisotropic tensor-product smooth uses this form: with one
parameter per margin, the fit can be rough in one direction and smooth
in another, whereas a single parameter smooths every direction alike.

## Usage

``` r
additive_penalty(
  matrices,
  map = NULL,
  blocks = 1L,
  readable = c("lambda", "sd"),
  link_lambda = linkfunctions7::log_link(),
  tol = 1e-10
)
```

## Arguments

- matrices:

  A list of at least one symmetric positive semidefinite matrix, all
  with the same number of rows. Each component is symmetrized. A
  component is rejected with an error that gives its position when it
  has an entry that is not finite (`NA`, `NaN` or infinite), when it
  differs from its transpose by more than `1e-8` times the larger of 1
  and its largest absolute entry, when its largest eigenvalue is not
  positive (a zero or negative semidefinite matrix), or when its
  smallest eigenvalue is below `-tol` times its largest. Under
  `blocks > 1` they are the components of one block.

- map:

  The matrix \\D\\, with as many rows as each component has and one
  column per coefficient, or `NULL` (the default) for the identity. A
  map is absorbed at construction, and each component becomes \\D^\top
  P_kD\\, so the stored `map` is always `NULL` and `n_coef` is the
  number of columns of the map. A map with an entry that is not finite
  is rejected, and so is a map that carries a component to the zero
  matrix, with an error that gives the position of that component.

- blocks:

  How many times each component is repeated blockwise, a whole number of
  at least 1. `1L`, the default, gives the ordinary penalty. A value
  above 1 cannot be combined with `map`.

- readable:

  `"lambda"` (the default) or `"sd"`, the quantities that
  [`penalty_readable()`](https://statmodels7.github.io/penalties7/reference/penalty_readable.md)
  reports. See the section on the readable quantities.

- link_lambda:

  The linkfunctions7 link that carries each smoothing parameter onto the
  whole real line,
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  by default. The same link is used for every parameter.

- tol:

  The relative eigenvalue tolerance, `1e-10` by default. It is used to
  reject a component whose smallest eigenvalue is below `-tol` times its
  largest, and to count the rank of the sum of the components after each
  has been divided by its largest absolute entry. The object stores
  `tol`, and
  [`additive_sum()`](https://statmodels7.github.io/penalties7/reference/additive_sum.md)
  passes it to
  [`additive_stable()`](https://statmodels7.github.io/penalties7/reference/additive_stable.md)
  as `r_tol`, so the rank of each group of components in the similarity
  transformation is counted at the same tolerance as the rank of the
  sum.

## Value

An
[`AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/AdditivePenalty.md)
object with one hyperparameter per component, named `lambda1`,
`lambda2`, ... The hyperparameters are bounded on \\(0, \infty)\\.

## The value

\$\$\rho(\beta;\lambda) = \tfrac{1}{2}\beta^\top S(\lambda)\beta -
\tfrac{1}{2}\log\mathrm{pdet}\\S(\lambda) + \tfrac{r}{2}\log 2\pi,
\qquad S(\lambda) = \sum\_{k} \lambda_k P_k,\$\$

with \\r\\ the rank of the sum. The hyperparameters are named `lambda1`,
`lambda2` and so on, in the order in which the components were given.
Each is positive and is carried onto the real line by `link_lambda`, a
log link by default.

## The derivatives of the log pseudo-determinant

With \\S^{+}\\ the pseudo-inverse of \\S(\lambda)\\, the derivatives
that a marginal criterion reads are traces against it:

\$\$\frac{\partial}{\partial\lambda_k}\log\mathrm{pdet}\\S =
\operatorname{tr}(S^{+}P_k), \qquad
\frac{\partial^{2}}{\partial\lambda_k\partial\lambda_l}
\log\mathrm{pdet}\\S = -\operatorname{tr}(S^{+}P_kS^{+}P_l).\$\$

Both formulas are exact.
[`additive_sum()`](https://statmodels7.github.io/penalties7/reference/additive_sum.md)
evaluates the traces, together with the log pseudo-determinant, by the
similarity transformation of
[`additive_stable()`](https://statmodels7.github.io/penalties7/reference/additive_stable.md),
and never forms \\S^{+}\\ itself.

## The rank of the sum

The rank and the null basis are computed once, at construction, from the
sum of the components after each has been divided by its largest
absolute entry, and the rank is the number of eigenvalues of that sum
above `tol` times the largest. Both are then fixed: the null space of
\\S(\lambda)\\ is the intersection of the null spaces of the components
at every positive \\\lambda\\. An eigenvalue count of the assembled
\\S(\lambda)\\ itself is correct only while the parameters are
comparable. Once they differ by many orders of magnitude, the smaller
contributions fall below any fixed relative tolerance and the count
drops, although the stored null basis still annihilates the sum to
rounding error.

## Parameters that differ by many orders of magnitude

The stored rank fixes how many eigenvalues enter the log
pseudo-determinant, but it does not make those eigenvalues computable
from the assembled sum. When the parameters differ by enough orders of
magnitude, the condition number of \\S(\lambda)\\ exceeds \\1/\epsilon\\
and the smallest eigenvalue on the range falls below the absolute
accuracy of one symmetric eigendecomposition, which is of order
\\d\\\epsilon\lVert S\rVert\_{2}\\. The log pseudo-determinant is
therefore not taken, in general, from one decomposition of the assembled
sum.
[`additive_stable()`](https://statmodels7.github.io/penalties7/reference/additive_stable.md),
which
[`additive_sum()`](https://statmodels7.github.io/penalties7/reference/additive_sum.md)
calls, partitions the components by the size of their contributions and
transforms by similarity, so that each decomposition resolves only the
spread within one group of components. The sum is decomposed as a whole
only when every contribution is at least
[`additive_tol()`](https://statmodels7.github.io/penalties7/reference/additive_tol.md)
times the largest, or when the dominant group by itself has the rank of
the sum. The identities of the transformation are exact, and the
accuracy of the value and of both derivative blocks is limited by the
conditioning of the dominant group at each level of the partition rather
than by the spread of the parameters. Where the weighted sum or one of
its reductions overflows, or where the rank identity that
[`additive_sum()`](https://statmodels7.github.io/penalties7/reference/additive_sum.md)
checks at every call is violated or cannot be evaluated, the log
pseudo-determinant and its derivatives are `NaN`, without an error. Such
points become more frequent as the spread and the number of components
grow.

## The quantities of a marginal criterion

[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
returns `TRUE`, because a sum of quadratic forms is a quadratic form, so
[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md),
[`penalty_rank()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md),
[`penalty_null_basis()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
and
[`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
all have methods for this class, and
[`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
runs its quadratic rows here (the row of the null basis only when the
null basis has columns). The matrix moves with one hyperparameter per
component, whereas the plain quadratic branch has a single scale, so the
log pseudo-determinant is not linear in the logarithm of any one
parameter. The row of
[`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
that checks it therefore compares its gradient with `numDeriv`, where
the plain quadratic branch checks a slope.

[`has_prox()`](https://statmodels7.github.io/penalties7/reference/has_prox.md)
returns `FALSE`. That function first tests whether the class registers a
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
method of its own, before it calls
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md),
and this class registers none.

`blocks = m` builds the penalty whose components are \\I_m \otimes
P_k\\, which is the form of a smooth repeated once per group. The sum is
\\I_m \otimes S_1(\lambda)\\, with \\S_1\\ the sum over one block, so
its rank is \\m\\ times the rank of the block, its log
pseudo-determinant is \\m\log\mathrm{pdet}\\S_1(\lambda)\\, and both
derivative blocks are \\m\\ times those of the block. The decomposition
is taken on one block, and the repeated components are stored as sparse
matrices, which is also the storage that
[`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
and
[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
return. `blocks` cannot be combined with `map`.

## The readable quantities

Setting `readable = "sd"` declares the smoothing parameters to be
precisions of a prior, and
[`penalty_readable()`](https://statmodels7.github.io/penalties7/reference/penalty_readable.md)
then reports the standard deviations \\\sigma_k = \lambda_k^{-1/2}\\,
named `sigma1`, `sigma2` and so on, in place of the parameters. The
penalty itself is unchanged.

## See also

[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
for one matrix with one scale,
[`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md)
for a matrix whose entries move with the hyperparameters,
[`penalty_logpdet.AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.AdditivePenalty.md)
for the log pseudo-determinant and its derivatives,
[`additive_sum()`](https://statmodels7.github.io/penalties7/reference/additive_sum.md)
for the computation that they share.

## Examples

``` r
# This is an anisotropic tensor smooth over a 4 x 4 grid, and the curvature
# along each margin is penalized separately.
P <- crossprod(diff(diag(4), differences = 2))
pen <- additive_penalty(list(kronecker(diag(4), P), kronecker(P, diag(4))))
pen@params
#> [1] "lambda1" "lambda2"
c(n_coef = pen@n_coef, rank = penalty_rank(pen))
#> n_coef   rank 
#>     16     12 

# Raising one parameter penalizes roughness along that margin alone.
set.seed(2)
b <- rnorm(16)
penalty_value(pen, b, list(lambda1 = 1, lambda2 = 1))
#> [1] 66.29988
penalty_value(pen, b, list(lambda1 = 1, lambda2 = 100))
#> [1] 2292.209

# With one component the value is that of the plain quadratic penalty.
P2 <- crossprod(diff(diag(5), differences = 2))
bb <- c(0.4, -1.1, 0.7, 0.2, -0.3)
penalty_value(additive_penalty(list(P2)), bb, list(lambda1 = 3)) -
  penalty_value(quadratic_penalty(P2), bb, list(lambda = 3))
#> [1] 3.552714e-15

# The stored rank stays fixed as the parameters spread apart, while an
# eigenvalue count of the assembled sum drops.
a <- additive_penalty(list(crossprod(diff(diag(5))), P2))
penalty_rank(a)
#> [1] 4
sapply(c(1, 1e10), function(r) {
  S <- penalty_matrix(a, list(lambda1 = 1, lambda2 = r))
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  sum(ev > 1e-10 * max(ev))
})
#> [1] 4 3

# blocks = m repeats each component without decomposing the repetition.
rep4 <- additive_penalty(list(crossprod(diff(diag(5))), P2), blocks = 4)
c(n_coef = rep4@n_coef, rank = penalty_rank(rep4))
#> n_coef   rank 
#>     20     16 
th <- list(lambda1 = 2, lambda2 = 0.5)
penalty_logpdet(rep4, th)$value - 4 * penalty_logpdet(a, th)$value
#> [1] 0

# readable = "sd" reports the standard deviations of the prior.
sdp <- additive_penalty(list(diag(2), diag(c(1, 2))), readable = "sd")
penalty_readable(sdp, list(lambda1 = 4, lambda2 = 0.25))$value
#> sigma1 sigma2 
#>    0.5    2.0 
```
