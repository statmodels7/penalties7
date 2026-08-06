# Construct a Quadratic Penalty

The exact negative log-density of the (possibly degenerate) Gaussian
prior \\\beta \sim N(0, (\lambda D'PD)^{-})\\: \$\$\rho =
\tfrac{\lambda}{2}\\(D\beta)'P(D\beta) - \tfrac{r}{2}\log\lambda +
\tfrac{r}{2}\log 2\pi - \tfrac{1}{2}\log\mathrm{pdet}(P),\$\$ with \\r\\
the rank of \\P\\. The constant is kept deliberately: dropping it, as
penalized-likelihood software does, makes joint hyperparameter
estimation degenerate.

## Usage

``` r
quadratic_penalty(
  P,
  map = NULL,
  link_lambda = linkfunctions7::log_link(),
  tol = 1e-10
)
```

## Arguments

- P:

  A symmetric positive semidefinite matrix, for instance a basis7 Gram
  matrix or a difference penalty \\D_k'D_k\\.

- map:

  The matrix \\D\\, or `NULL` (default) for the identity.

- link_lambda:

  The link carrying \\\lambda\\; defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

- tol:

  The relative eigenvalue tolerance of the rank rule.

## Value

An object of class `QuadraticPenalty`.

## Details

One eigendecomposition at construction fixes the rank by the relative
rule \\\mathrm{ev} \> \mathrm{tol}\cdot\max(\mathrm{ev})\\ – a statement
about the matrix rather than about whichever arithmetic is later
performed on it – and stores the exact null basis of \\D'PD\\, so that
membership questions never go through a rank recomputed from an
assembled sum. Every quantity is then closed form in \\\lambda\\.

## Examples

``` r
# a second-difference penalty on five coefficients: rank 3, an improper
# prior whose null space is the linear polynomials
P <- crossprod(diff(diag(5), differences = 2))
pen <- quadratic_penalty(P)
penalty_rank(pen)
#> [1] 3
is_proper(pen)
#> [1] FALSE
penalty_value(pen, rnorm(5), list(lambda = 2))
#> [1] 9.451446
```
