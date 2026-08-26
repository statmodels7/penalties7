# Value of a Penalty

Returns the scalar \\\rho(D\beta; \theta)\\, the amount a set of
coefficients is penalized by. The normalizing constant is included
whenever the penalty is proper, so the value is exactly the negative
log-density of the prior and can be added to a negative log-likelihood
without a further term. An improper penalty returns the bare \\\rho\\;
ask
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
which you have.

## Usage

``` r
penalty_value(pen, beta, theta, ...)
```

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object of any branch.

- beta:

  A numeric vector of length `pen@n_coef`. Coerced with
  [`as.numeric()`](https://rdrr.io/r/base/numeric.html) in the generic,
  so an integer vector or a one-column matrix is accepted.

- theta:

  A named list of hyperparameter values, or a named numeric vector
  carrying the same, holding every name in `pen@params`. Reordered and
  checked against the penalty's open bounds before dispatch. Pass
  [`list()`](https://rdrr.io/r/base/list.html) for a penalty with no
  hyperparameters.

- ...:

  Passed to methods. No shipped method reads it.

## Value

A single number. Never `NA` for an argument inside the bounds; a penalty
whose value is infinite at some \\\beta\\ returns `Inf`.

## Why the constant is kept

Dropping the normalizing constant makes no difference to a fit at a
fixed hyperparameter, since it does not depend on \\\beta\\. It makes
the hyperparameter itself unestimable: with the constant dropped,
driving \\\lambda\\ to zero drives the penalty to zero and the joint
maximum runs away. With it kept, a proper penalty is a density in
\\\beta\\ for every \\\lambda\\, so the joint objective has an interior
maximum and a marginal criterion has something to expand around.
Estimating the degrees of freedom of a heavy-tailed prior needs the
constant for the same reason.

## The map

\\\rho\\ is evaluated at \\t = D\beta\\, not at \\\beta\\. With the map
`NULL`, the default of every constructor, \\t = \beta\\ and no
arithmetic is done. A map of \\m\\ rows means \\\rho\\ sees \\m\\
numbers, so a second-difference map penalizes curvature and leaves the
level and the slope alone.

## See also

[`penalty_gradient()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
and
[`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
for the coefficient derivatives,
[`penalty_grad_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
for the hyperparameter ones,
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
for whether the constant is there,
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
for where the value is not differentiable.

## Examples

``` r
# A ridge over three coefficients, with the Gaussian constant kept.
pen <- quadratic_penalty(diag(3))
penalty_value(pen, c(1, 0, -1), list(lambda = 2))
#> [1] 3.717095

# Which is exactly the negative log-density of a N(0, 1/lambda) prior.
-sum(stats::dnorm(c(1, 0, -1), sd = 1 / sqrt(2), log = TRUE))
#> [1] 3.717095

# The constant is why the value moves with lambda even at beta = 0.
sapply(c(0.5, 1, 2, 8),
       function(l) penalty_value(pen, c(0, 0, 0), list(lambda = l)))
#> [1]  3.7965364  2.7568156  1.7170948 -0.3623467

# A map sends rho the second differences, so only curvature is charged for.
# A straight line contributes nothing to the quadratic part and the value
# falls to the constant, the same value the zero vector gives.
curve <- quadratic_penalty(diag(2), map = diff(diag(4), differences = 2))
penalty_value(curve, c(1, 2, 3, 4), list(lambda = 5))
#> [1] 0.2284392
penalty_value(curve, c(0, 0, 0, 0), list(lambda = 5))
#> [1] 0.2284392
penalty_value(curve, c(1, 2, 4, 8), list(lambda = 5))
#> [1] 12.72844
```
