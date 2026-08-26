# Value of an Additive Penalty

Returns the negative log-density of the Gaussian prior whose precision
is the weighted sum \\S(\lambda) = \sum_k \lambda_k P_k\\, with the
normalizing constant taken over the range of the sum.

## Arguments

- pen:

  An
  [`AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/AdditivePenalty.md)
  object.

- beta:

  A numeric vector of length `pen@n_coef`, already coerced by the
  generic.

- theta:

  A named list of `lambda1`, `lambda2`, ..., already aligned and
  bound-checked by the generic.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A single number.

## Details

\$\$\rho(\beta;\lambda) = \tfrac{1}{2}\beta^\top S(\lambda)\beta -
\tfrac{1}{2}\log\mathrm{pdet}\\S(\lambda) + \tfrac{r}{2}\log 2\pi,\$\$

with \\r\\ the stored rank. One eigendecomposition of the sum, through
[`additive_sum()`](https://statmodels7.github.io/penalties7/reference/additive_sum.md),
supplies both the quadratic form's matrix and the log
pseudo-determinant.

With a single component the value is the plain quadratic penalty's, and
the two agree to `3.6e-15`.

## See also

[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
for the construction,
[`penalty_gradient.AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.AdditivePenalty.md)
for the coefficient derivatives,
[`additive_sum()`](https://statmodels7.github.io/penalties7/reference/additive_sum.md)
for the decomposition it reads.

## Examples

``` r
P1 <- crossprod(diff(diag(5)))
P2 <- crossprod(diff(diag(5), differences = 2))
pen <- additive_penalty(list(P1, P2))
b <- c(0.4, -1.1, 0.7, 0.2, -0.3)

penalty_value(pen, b, list(lambda1 = 2, lambda2 = 0.5))
#> [1] 10.94012

# A constant vector lies in the intersected null space, so only the
# normalizing constant remains.
penalty_value(pen, rep(1, 5), list(lambda1 = 2, lambda2 = 0.5))
#> [1] 0.9051224
penalty_value(pen, rep(0, 5), list(lambda1 = 2, lambda2 = 0.5))
#> [1] 0.9051224
```
