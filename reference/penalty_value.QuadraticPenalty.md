# Value of a Quadratic Penalty

Returns the negative log-density of the Gaussian prior with precision
\\\lambda D'PD\\, evaluated by adding the constant to the quadratic
form. Every term is closed in \\\lambda\\, the rank and the log
pseudo-determinant of \\P\\ having been fixed at construction.

## Arguments

- pen:

  A
  [`QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/QuadraticPenalty.md)
  object.

- beta:

  A numeric vector of length `pen@n_coef`, already coerced by the
  generic.

- theta:

  A list holding `lambda`, already aligned and bound-checked by the
  generic.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A single number.

## Details

\$\$\rho(\beta; \lambda) = \frac{\lambda}{2}\\(D\beta)'P(D\beta) -
\frac{r}{2}\log\lambda + \frac{r}{2}\log 2\pi -
\frac{1}{2}\log\mathrm{pdet}(P),\$\$

with \\r\\ the stored rank. At \\\beta = 0\\ the value is the constant
alone and still moves with \\\lambda\\, so the smoothing parameter stays
estimable from a joint objective.

## See also

[`penalty_gradient.QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.QuadraticPenalty.md)
for the coefficient derivatives,
[`penalty_logpdet.QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.QuadraticPenalty.md)
for the constant's own piece,
[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
for the constructor.

## Examples

``` r
# A full-rank penalty is a proper prior, so the value is a log-density.
pen <- quadratic_penalty(diag(3))
b <- c(0.4, -1.1, 0.7)
all.equal(penalty_value(pen, b, list(lambda = 2)),
          -sum(stats::dnorm(b, sd = 1 / sqrt(2), log = TRUE)))
#> [1] TRUE

# The constant alone, at the origin, and how it moves with lambda.
sapply(c(0.5, 2, 8),
       function(l) penalty_value(pen, c(0, 0, 0), list(lambda = l)))
#> [1]  3.7965364  1.7170948 -0.3623467
```
