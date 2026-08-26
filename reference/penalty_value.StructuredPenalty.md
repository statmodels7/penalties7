# Value of a Structured Penalty

Returns the negative log-density of the Gaussian prior the structure
describes, evaluated in the precision whichever of the two matrices the
structure supplies.

## Arguments

- pen:

  A
  [`StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/StructuredPenalty.md)
  object.

- beta:

  A numeric vector of length `pen@n_coef`, already coerced by the
  generic.

- theta:

  A named list of the structure's free values, already aligned by the
  generic. The bounds are the whole line, so nothing is rejected here.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A single number.

## Details

\$\$\rho(\beta; \theta) = \tfrac{1}{2}\\\beta'\Omega(\theta)\beta -
\tfrac{1}{2}\log\mathrm{pdet}\\\Omega(\theta) + \tfrac{r}{2}\log
2\pi,\$\$

with \\r\\ the structure's rank, and the log-determinant negated where
the structure is a covariance. For a full-rank structure this is exactly
a multivariate normal log-density with mean zero.

## See also

[`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md)
for the construction,
[`penalty_gradient.StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.StructuredPenalty.md)
for the coefficient derivatives,
[`struct_omega()`](https://statmodels7.github.io/penalties7/reference/struct_omega.md)
for the transport.

## Examples

``` r
# A covariance structure gives exactly a multivariate normal log-density.
s <- parameters7::ar1(4, role = "covariance")
pen <- structured_penalty(s)
b <- c(0.3, -0.1, 0.4, 0.2)
S <- unclass(parameters7::param_value(s, c(0.2, 0.5)))
penalty_value(pen, b, list(log_scale = 0.2, z_rho = 0.5)) -
  (0.5 * sum(b * solve(S, b)) +
     0.5 * as.numeric(determinant(S, logarithm = TRUE)$modulus) +
     2 * log(2 * pi))
#> [1] 0
```
