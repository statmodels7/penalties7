# Value of a Separable Penalty

Returns \\-\sum_i \log f(b_i;\theta)\\, the negative log-density of the
parent summed over the blocks of \\D\beta\\. Because the parent supplies
its own normalizing constant, the value is exactly the negative log
prior density and needs nothing added.

## Arguments

- pen:

  A
  [`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
  object.

- beta:

  A numeric vector of length `pen@n_coef`, already coerced by the
  generic.

- theta:

  A named list of the parent's free parameters, already aligned and
  bound-checked by the generic against the parent's own bounds.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A single number.

## Details

The map is applied, the result reshaped into the parent's argument by
[`dp_arg()`](https://statmodels7.github.io/penalties7/reference/dp_arg.md),
and one call to
[`distributions7::distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.html)
with `log = TRUE` gives every block's contribution at once. A univariate
parent sees the vector; a \\p\\-variate one sees a matrix of one row per
block.

## See also

[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md)
for the construction,
[`penalty_gradient.DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.DistribPenalty.md)
for the coefficient derivatives,
[`ridge_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md)
and its siblings for the named instances.

## Examples

``` r
# A Gaussian prior at zero, which is the negative log-density of one.
pen <- distrib_penalty(
  distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0),
  n_coef = 3)
b <- c(1, 0, -1)
penalty_value(pen, b, list(sigma = 2))
#> [1] 5.086257
-sum(stats::dnorm(b, sd = 2, log = TRUE))
#> [1] 5.086257

# The lasso's value is lambda times the L1 norm plus -q log(lambda/2).
penalty_value(lasso_penalty(n_coef = 3), b, list(lambda = 0.5))
#> [1] 5.158883
0.5 * sum(abs(b)) - 3 * log(0.5 / 2)
#> [1] 5.158883
```
