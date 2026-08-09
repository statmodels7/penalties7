# Is a Penalty Quadratic?

`TRUE` only for
[`quadratic_penalty`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md),
whose matrix, rank, null basis and log pseudo-determinant the
marginal-likelihood generics expose.

## Usage

``` r
is_quadratic(pen, ...)
```

## Arguments

- pen:

  A
  [`penalty`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- ...:

  Passed to methods.

## Value

A single logical.

## Details

A quadratic penalty carries a fixed matrix \\P\\ and one smoothing
parameter \\\lambda\\, and is the negative log-density of the improper
Gaussian prior with precision \\\lambda P\\ on \\D\beta\\:

\$\$\rho(\beta; \lambda) = \tfrac{\lambda}{2} (D\beta)^\top P (D\beta) -
\tfrac{1}{2}\log^{+}\lvert \lambda P \rvert + \tfrac{r}{2}\log(2\pi),
\qquad \log^{+}\lvert \lambda P \rvert = r \log \lambda + \log^{+}\lvert
P \rvert,\$\$

with \\\log^{+}\\ the log pseudo-determinant and \\r =
\operatorname{rank}(P)\\. Those are the quantities
[`penalty_matrix`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md),
[`penalty_rank`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md),
[`penalty_null_basis`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
and
[`penalty_logpdet`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
report and a REML or marginal-likelihood criterion needs; a penalty for
which this is `FALSE` has no such matrix and those generics reject.

## Examples

``` r
is_quadratic(quadratic_penalty(diag(2)))
#> [1] TRUE
is_quadratic(ridge_penalty())
#> [1] FALSE
```
