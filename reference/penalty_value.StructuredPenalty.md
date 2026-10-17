# Value of a Structured Penalty

Returns the negative log-density of the Gaussian prior whose precision
\\\Omega(\theta)\\ is the matrix of the structure. A structure built
with
[`parameters7::inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.html)
supplies the inverse of the matrix that it wraps, so that matrix is then
the covariance.

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
  generic. The bounds are the whole line, so this method rejects no
  value.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A single number.

## Details

\$\$\rho(\beta; \theta) = \tfrac{1}{2}\\\beta'\Omega(\theta)\beta -
\tfrac{1}{2}\log\mathrm{pdet}\\\Omega(\theta) + \tfrac{r}{2}\log
2\pi,\$\$

with \\r\\ the structure's rank. For a full-rank structure this is the
negative log-density of a multivariate normal distribution with mean
zero and precision \\\Omega(\theta)\\.

## See also

[`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md)
for the construction,
[`penalty_gradient.StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.StructuredPenalty.md)
for the coefficient derivatives,
[`struct_omega()`](https://statmodels7.github.io/penalties7/reference/struct_omega.md)
for the helpers that read the precision from the structure.

## Examples

``` r
# A structure wrapped in inverse_of() places the AR(1) matrix S on the
# covariance side, and the value is the negative log-density of N(0, S).
s <- parameters7::inverse_of(parameters7::ar1(4))
pen <- structured_penalty(s)
b <- c(0.3, -0.1, 0.4, 0.2)
S <- unclass(parameters7::param_value(parameters7::ar1(4), c(0.2, 0.5)))
penalty_value(pen, b, list(log_scale = 0.2, z_rho = 0.5)) -
  (0.5 * sum(b * solve(S, b)) +
     0.5 * as.numeric(determinant(S, logarithm = TRUE)$modulus) +
     2 * log(2 * pi))
#> [1] 0
```
