# Construct a Sum of Quadratic Penalties

The gaussian prior whose precision is a weighted sum of fixed matrices,
\$\$\rho(\beta;\lambda) = \tfrac{1}{2}\beta^\top S(\lambda)\beta -
\tfrac{1}{2}\log\mathrm{pdet}\\S(\lambda) + \tfrac{r}{2}\log 2\pi,
\qquad S(\lambda) = \sum\_{k} \lambda_k P_k,\$\$ with one smoothing
parameter per component. This is what a tensor-product smooth needs: one
parameter per margin, so that the fit may be rough in one direction and
smooth in another.

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

  A list of symmetric positive semidefinite matrices of the same
  dimension.

- map:

  The matrix \\D\\, or `NULL` (default) for the identity.

- link_lambda:

  The link carrying each smoothing parameter to the unconstrained scale.
  Defaults to the log.

- tol:

  The relative tolerance below which an eigenvalue counts as zero.

## Value

An object of class
[`AdditivePenalty`](https://statmodels7.github.io/penalties7/reference/AdditivePenalty.md).

## Details

[`quadratic_penalty`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
carries a single matrix and one scale, which forces every direction to
be smoothed alike. Here the components keep their own parameters, and
the quantities a marginal criterion reads follow from one
eigendecomposition of the sum:
\$\$\frac{\partial}{\partial\lambda_k}\log\mathrm{pdet}\\S =
\operatorname{tr}(S^{+}P_k), \qquad
\frac{\partial^{2}}{\partial\lambda_k\partial\lambda_l}
\log\mathrm{pdet}\\S = -\operatorname{tr}(S^{+}P_kS^{+}P_l).\$\$

**The rank is not read off the sum.** Counting the eigenvalues of
\\S(\lambda)\\ above a tolerance gives the right answer only while the
parameters are comparable: once they differ by orders of magnitude the
small contributions sink below any fixed tolerance and are counted as
zeros, so the rank appears to fall as the fit is smoothed. The null
space of a sum of positive semidefinite matrices is the intersection of
their null spaces, which does not depend on the parameters at all, so
the rank is fixed once at construction from the components stacked and
individually normalized.

## See also

[`quadratic_penalty`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)

## Examples

``` r
# curvature in two directions, penalized separately
P1 <- crossprod(diff(diag(4), differences = 2))
pen <- additive_penalty(list(kronecker(diag(4), P1),
                             kronecker(P1, diag(4))))
pen@params
#> [1] "lambda1" "lambda2"
penalty_value(pen, rnorm(16), list(lambda1 = 1, lambda2 = 100))
#> [1] 1895.369
```
