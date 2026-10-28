# Value of a Separable Penalty

Returns \\-\sum_i \log f(b_i;\theta)\\, the negative log-density of the
parent summed over the blocks \\b_i\\ of \\D\beta\\. Because the parent
supplies its own normalizing constant, the value is the negative log
prior density of \\D\beta\\. Under a square nonsingular map, diagonal or
not, \\\log\lvert\det D\rvert\\ is subtracted, which makes the value the
negative log-density of \\\beta\\. Under a map that is not square, or is
singular, the function does not subtract a constant.

## Arguments

- penalty:

  A
  [`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
  object.

- beta:

  A numeric vector of length `penalty@n_coef`, already coerced by the
  generic.

- theta:

  A named list of the parent's free parameters, already aligned and
  bound-checked by the generic against the parent's own bounds.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A single number.

## Details

The map is applied, the result is reshaped into the argument of the
parent by
[`dp_arg()`](https://statmodels7.github.io/penalties7/reference/dp_arg.md),
and one call to
[`distributions7::distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.html)
with `log = TRUE` gives the contribution of every block at once. A
univariate parent receives the vector, and a \\p\\-variate one a matrix
of one row per block.

The log-determinant is read by
[dp_log_jacobian()](https://statmodels7.github.io/penalties7/reference/dp_diag_entries.md):
it is computed from the entries of a diagonal map, and for any other map
it is the `log_jacobian` stored when the penalty was built.

## See also

[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md)
for the construction,
[`penalty_gradient.DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.DistribPenalty.md)
for the coefficient derivatives,
[`ridge_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md)
and its siblings for the named instances.

## Examples

``` r
# A Gaussian prior at zero: the value is its negative log-density.
g0 <- distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0)
pen <- distrib_penalty(g0, n_coef = 3)
b <- c(1, 0, -1)
penalty_value(pen, b, list(sigma = 2))
#> [1] 5.086257
-sum(stats::dnorm(b, sd = 2, log = TRUE))
#> [1] 5.086257

# The value of the lasso is lambda times the L1 norm minus q log(lambda/2).
penalty_value(lasso_penalty(n_coef = 3), b, list(lambda = 0.5))
#> [1] 5.158883
0.5 * sum(abs(b)) - 3 * log(0.5 / 2)
#> [1] 5.158883

# A diagonal map subtracts the sum of the log-moduli of its entries, and a
# map with fewer rows than columns subtracts nothing.
dg <- distrib_penalty(g0, map = Matrix::Diagonal(x = c(2, 0.5, 4)))
penalty_value(dg, b, list(sigma = 2))
#> [1] 5.949963
-sum(stats::dnorm(c(2, 0, -4), sd = 2, log = TRUE)) - log(2 * 0.5 * 4)
#> [1] 5.949963
sh <- distrib_penalty(g0, map = diff(diag(3)))
penalty_value(sh, b, list(sigma = 2))
#> [1] 3.474171
-sum(stats::dnorm(diff(b), sd = 2, log = TRUE))
#> [1] 3.474171
```
