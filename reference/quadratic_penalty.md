# Construct a Quadratic Penalty

Builds the penalty \\\tfrac{\lambda}{2}(D\beta)'P(D\beta)\\ together
with the constant that makes it a density, so the value returned is
exactly the negative log-density of the Gaussian prior with precision
\\\lambda D'PD\\. This is the penalty of a spline smooth, a ridge and a
Gaussian random effect: one fixed matrix that weights the directions of
the coefficients, and one smoothing parameter that scales the weights.

## Usage

``` r
quadratic_penalty(
  P,
  map = NULL,
  blocks = 1L,
  link_lambda = linkfunctions7::log_link(),
  tol = 1e-10
)
```

## Arguments

- P:

  A symmetric positive semidefinite matrix, for instance a basis7 Gram
  matrix or a difference penalty \\D_k'D_k\\. A matrix with an entry
  that is not finite (`NA`, `NaN` or infinite) is rejected. Symmetry is
  checked to a tolerance of `1e-8` times the larger of 1 and the largest
  absolute entry of `P`, and the matrix is symmetrized before use. A
  matrix whose largest eigenvalue is not positive (a zero or negative
  semidefinite matrix) is rejected, and so is one with an eigenvalue
  below `-tol` times the largest.

- map:

  The matrix \\D\\, with as many rows as \\P\\ and one column per
  coefficient, or `NULL` (the default) for the identity, and then
  `n_coef` is `nrow(P)`. A Matrix object keeps its own storage. A map
  with an entry that is not finite is rejected, and so is a map for
  which \\D'PD\\ has no positive eigenvalue; a map that annihilates
  \\P\\ only up to rounding is not detected (see the details).

- blocks:

  How many times \\P\\ is repeated blockwise, a whole number of at
  least 1. `1L`, the default, is the ordinary single-block penalty.
  Above 1 it must be given without a `map`.

- link_lambda:

  The linkfunctions7 link carrying \\\lambda\\ onto the whole real line.
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  by default, \\\lambda\\ being positive.

- tol:

  The relative eigenvalue tolerance of the rank rule, `1e-10` by
  default. An eigenvalue of \\D'PD\\ at or below `tol` times the largest
  counts as zero, so its direction joins the null space and is not
  penalized. The same tolerance bounds the negative eigenvalues of \\P\\
  that are accepted as rounding errors. The object stores `tol`, and
  [`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
  counts the rank at the same tolerance.

## Value

A
[`QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/QuadraticPenalty.md)
object with one hyperparameter, `lambda`, bounded on \\(0, \infty)\\.

## The value

With \\r\\ the rank of \\D'PD\\,

\$\$\rho(\beta; \lambda) = \tfrac{\lambda}{2}\\(D\beta)'P(D\beta) -
\tfrac{r}{2}\log\lambda + \tfrac{r}{2}\log 2\pi -
\tfrac{1}{2}\log\mathrm{pdet}(D'PD).\$\$

The last three terms are the normalizing constant of the prior, taken
over the range of \\D'PD\\ alone when that matrix is rank deficient. The
constant is the one of the matrix applied to \\\beta\\, so it carries
the map. For a diagonal map with non-zero entries \\d_j\\ and a
full-rank \\P\\, \\\log\mathrm{pdet}(D'PD) = \log\det P +
2\sum_j\log\lvert d_j\rvert\\, and \\\sum_j\log\lvert d_j\rvert\\ is the
log Jacobian of \\\beta \mapsto D\beta\\. The value is then the negative
log-density of a proper Gaussian prior on \\\beta\\, and a marginal
criterion built on it does not change when a covariate is multiplied by
a factor and the map multiplies the coefficient of that covariate by the
same factor. The constant also makes \\\lambda\\ estimable from a joint
objective: without the term \\-\tfrac{r}{2}\log\lambda\\ the penalty
decreases to zero with \\\lambda\\, and the joint maximum moves to the
boundary \\\lambda \to 0\\.

## The rank and the null basis

The rank is counted at construction by the relative rule \\\mathrm{ev}
\> \mathtt{tol} \cdot \max(\mathrm{ev})\\ on the eigenvalues of
\\D'PD\\, and the same mask selects the orthonormal null basis among the
eigenvectors, so the null basis always has \\q - r\\ columns, with \\q\\
the number of coefficients. Without a map the eigenvalues come from a
values-only decomposition of \\P\\ and the eigenvectors from a second
decomposition of the same matrix; both are sorted in decreasing order,
so the mask selects the matching columns. With a map the rank, the log
pseudo-determinant and the null basis all come from the decomposition of
\\D'PD\\.

Fixing the rank on each matrix keeps it stable when quadratic penalties
are combined. A rank recounted from an assembled sum of several
penalties falls as their smoothing parameters spread apart, while the
true null space is the intersection of the null spaces of the components
and does not depend on the smoothing parameters.

The rule is relative to the largest eigenvalue of \\D'PD\\ itself. A map
that annihilates \\P\\ only up to rounding, so that every entry of
\\D'PD\\ is a rounding error relative to \\P\\, is therefore not
detected: the largest of those rounding errors passes the rule, and the
penalty is built with a positive rank and a log pseudo-determinant made
of rounding errors. The constructor rejects the map only when no
eigenvalue of \\D'PD\\ is positive, as when \\D'PD\\ is exactly the zero
matrix.

## Blockwise repetition

`blocks = m` builds the penalty of \\I_m \otimes P\\ from \\P\\ alone,
which is the form that one copy of a smooth per level of a factor takes.
The eigenvalues of \\I_m \otimes P\\ are those of \\P\\ repeated \\m\\
times, so the rank is \\mr\\, the log pseudo-determinant is
\\m\log\mathrm{pdet}(P)\\ and the null space is \\I_m \otimes N\\, with
\\N\\ the null basis of \\P\\. Only \\P\\ is decomposed, and the large
matrix is stored in sparse form. The result has the same value,
gradient, Hessian, rank, null basis and log pseudo-determinant as the
penalty built from the assembled \\I_m \otimes P\\, up to rounding.

**Under `blocks > 1` the matrix-valued returns are `dgCMatrix`**, not
base matrices:
[`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md),
[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md),
[`penalty_null_basis()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
and the matrices of
[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
all keep the sparse storage, so the assembled matrix is never densified.
A consumer that writes these results into a block of its own information
matrix has to accept both classes.

`blocks` does not combine with `map`, and the constructor rejects the
pair: a map mixes the blocks, so \\D'(I_m \otimes P)D\\ is in general
not block diagonal and does not have the structure that `blocks` builds.

## See also

[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
for a sum of quadratics with a smoothing parameter each,
[`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md)
for a matrix that moves with several hyperparameters,
[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md)
for the separable branch,
[`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
for the constant that a marginal criterion reads,
[`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
to verify the result.

## Examples

``` r
# Second differences over five coefficients: rank 3, and a null space
# holding the straight lines, so this is an improper prior.
P <- crossprod(diff(diag(5), differences = 2))
pen <- quadratic_penalty(P)
penalty_rank(pen)
#> [1] 3
is_proper(pen)
#> [1] FALSE

# A straight line lies in the null space and has zero gradient.
penalty_gradient(pen, c(1, 2, 3, 4, 5), list(lambda = 2))
#> [1] 0 0 0 0 0

# A full-rank penalty is a proper Gaussian prior, and its value is the
# negative log-density of one.
ridge <- quadratic_penalty(diag(3))
b <- c(0.4, -1.1, 0.7)
all.equal(penalty_value(ridge, b, list(lambda = 2)),
          -sum(stats::dnorm(b, sd = 1 / sqrt(2), log = TRUE)))
#> [1] TRUE

# A map penalizes its image. Here only the second differences of a
# six-vector are penalized, through a 4 x 6 map.
D <- diff(diag(6), differences = 2)
curved <- quadratic_penalty(diag(4), map = D)
curved@n_coef
#> [1] 6
penalty_value(curved, 1:6, list(lambda = 1))
#> [1] 1.348774

# blocks = m repeats P without assembling I_m (x) P.
blocked <- quadratic_penalty(P, blocks = 4)
c(blocked@n_coef, penalty_rank(blocked))
#> [1] 20 12
all.equal(penalty_rank(blocked), 4L * penalty_rank(pen))
#> [1] TRUE
```
