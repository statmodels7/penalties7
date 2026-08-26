# Construct a Quadratic Penalty

Builds the penalty \\\tfrac{\lambda}{2}(D\beta)'P(D\beta)\\ together
with the constant that makes it a density, so the value returned is
exactly the negative log-density of the Gaussian prior with precision
\\\lambda D'PD\\. This is the penalty behind a spline smooth, a ridge
and a Gaussian random effect: one fixed matrix saying which directions
are expensive, and one smoothing parameter saying how expensive.

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
  matrix or a difference penalty \\D_k'D_k\\. Symmetry is checked to a
  relative tolerance of `1e-8` and the matrix is symmetrized before use;
  a zero matrix is rejected.

- map:

  The matrix \\D\\, with as many rows as \\P\\ and one column per
  coefficient, or `NULL` (the default) for the identity, and then
  `n_coef` is `nrow(P)`. A Matrix object keeps its own storage.

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
  default. An eigenvalue at or below `tol` times the largest counts as
  zero, so its direction joins the null space and is not penalized.

## Value

A
[`QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/QuadraticPenalty.md)
object with one hyperparameter, `lambda`, bounded on \\(0, \infty)\\.

## The value

With \\r\\ the rank of \\P\\,

\$\$\rho(\beta; \lambda) = \tfrac{\lambda}{2}\\(D\beta)'P(D\beta) -
\tfrac{r}{2}\log\lambda + \tfrac{r}{2}\log 2\pi -
\tfrac{1}{2}\log\mathrm{pdet}(P).\$\$

The last three terms are the normalizing constant of the prior, taken
over the range of \\P\\ alone when \\P\\ is deficient.
Penalized-likelihood software usually drops them, and dropping them
makes \\\lambda\\ unestimable: with no \\-\tfrac{r}{2}\log\lambda\\ the
penalty falls to zero as \\\lambda\\ does and the joint maximum runs
away.

## The rank, and why it is fixed at construction

One eigendecomposition fixes the rank by the relative rule \\\mathrm{ev}
\> \mathtt{tol} \cdot \max(\mathrm{ev})\\, and stores the exact null
basis of \\D'PD\\. That is a statement about the matrix, not about
whichever arithmetic is later performed on it, and it is what keeps the
answer stable: a rank recounted from an assembled sum of several
penalties falls as their smoothing parameters spread apart, while the
true null space is the intersection of the components' and does not
move.

## Blockwise repetition

`blocks = m` builds the penalty of \\I_m \otimes P\\ from \\P\\ alone,
the shape one copy of a smooth per level of a factor takes. The
eigenvalues of \\I_m \otimes P\\ are \\P\\'s repeated \\m\\ times, so
the rank is \\mr\\, the log pseudo-determinant is
\\m\log\mathrm{pdet}(P)\\ and the null space is \\I_m \otimes N\\: the
large matrix is never decomposed. Measured on second differences over
ten coefficients at \\m = 200\\, so an assembled matrix of \\2000 \times
2000\\:

|                    |               |                |
|--------------------|---------------|----------------|
|                    | assembled     | `blocks = 200` |
| eigendecomposition | 4.79 s        | 2.9e-05 s      |
| whole construction | 5.82 s        | 0.001 s        |
| stored matrix      | 30.5 MB dense | 0.11 MB sparse |

at a density of 0.0022. The two penalties agree exactly: the value to
1.4e-14, the gradient and the Hessian to 0, and the same rank, null
basis and log pseudo-determinant.

**Under `blocks > 1` the matrix-valued returns are `dgCMatrix`**, not
base matrices:
[`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md),
[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md),
[`penalty_null_basis()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
and
[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
all keep the sparse storage. A branch built to avoid the assembled form
would defeat itself by densifying at its own boundary. Consumers that
write these into a block of their own information need to accept both
classes.

`blocks` does not combine with `map`, and the constructor rejects the
pair: a map mixes the blocks, so \\D'(I_m \otimes P)D\\ is not block
diagonal with the structure `blocks` names.

## See also

[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
for a sum of quadratics with a smoothing parameter each,
[`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md)
for a matrix that moves with several hyperparameters,
[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md)
for the separable branch,
[`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
for the constant a marginal criterion reads,
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

# A map penalizes what it selects. Here only the second differences of a
# six-vector are charged for, through a 4 x 6 map.
D <- diff(diag(6), differences = 2)
curved <- quadratic_penalty(diag(4), map = D)
curved@n_coef
#> [1] 6
penalty_value(curved, 1:6, list(lambda = 1))
#> [1] 3.675754

# blocks = m repeats P without assembling I_m (x) P.
blocked <- quadratic_penalty(P, blocks = 4)
c(blocked@n_coef, penalty_rank(blocked))
#> [1] 20 12
all.equal(penalty_rank(blocked), 4L * penalty_rank(pen))
#> [1] TRUE
```
