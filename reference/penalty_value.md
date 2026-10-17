# Value of a Penalty

Returns the scalar \\\rho(D\beta; \theta)\\, the amount by which a set
of coefficients is penalized. The normalizing constant is included
whenever the penalty is proper, so the value is exactly the negative
log-density of the prior and can be added to a negative log-likelihood
without a further term.
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
reports whether this is the case. An improper penalty keeps the
constant, if any, that its branch defines. A rank-deficient quadratic,
additive or structured penalty keeps a log pseudo-determinant, which is
taken over the range of its matrix. A separable penalty under a map that
is not square, or is singular, returns \\-\sum_i \log f(t_i; \theta)\\
without a Jacobian term. SCAD and MCP return the bare \\\rho\\.

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

  A finite numeric vector of length `pen@n_coef`. The generic coerces it
  with [`as.numeric()`](https://rdrr.io/r/base/numeric.html), so an
  integer vector or a one-column matrix is accepted. A vector of another
  length signals an error.

- theta:

  A named list of hyperparameter values, or a named numeric vector
  carrying the same, holding every name in `pen@params`. The generic
  reorders it and checks it against the penalty's open bounds before
  dispatch. Pass [`list()`](https://rdrr.io/r/base/list.html) for a
  penalty with no hyperparameters.

- ...:

  Passed to methods. The shipped methods do not read it.

## Value

A single number, which is not `NA` when `beta` is finite and `theta`
lies inside the bounds. The value is defined only for a finite `beta`:
at an infinite coefficient a quadratic penalty returns `NaN`.

## The normalizing constant

Dropping the normalizing constant makes no difference to a fit at a
fixed hyperparameter, since it does not depend on \\\beta\\. It does
make the hyperparameter itself impossible to estimate: if the constant
is dropped, driving \\\lambda\\ to zero drives the penalty to zero and
the joint maximum runs away. If the constant is kept, a proper penalty
is a density in \\\beta\\ for every \\\lambda\\, so the joint objective
has an interior maximum, and a marginal criterion has a mode around
which it can be expanded. The degrees of freedom of a heavy-tailed prior
are estimable for the same reason.

## The map

\\\rho\\ is evaluated at \\t = D\beta\\, not at \\\beta\\. If the map is
`NULL`, which is the default of every constructor that takes a map, then
\\t = \beta\\ and the function does not apply the map. A map of \\m\\
rows gives \\\rho\\ an argument of length \\m\\, so a second-difference
map penalizes curvature and leaves the level and the slope unpenalized.

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
# A ridge over three coefficients includes the Gaussian constant.
pen <- quadratic_penalty(diag(3))
penalty_value(pen, c(1, 0, -1), list(lambda = 2))
#> [1] 3.717095

# This value is exactly the negative log-density of a N(0, 1/lambda) prior.
-sum(stats::dnorm(c(1, 0, -1), sd = 1 / sqrt(2), log = TRUE))
#> [1] 3.717095

# Because of the constant, the value moves with lambda even at beta = 0.
sapply(c(0.5, 1, 2, 8),
       function(l) penalty_value(pen, c(0, 0, 0), list(lambda = l)))
#> [1]  3.7965364  2.7568156  1.7170948 -0.3623467

# Under a second-difference map rho is evaluated at the second differences,
# so a straight line contributes nothing to the quadratic part and the
# value equals the constant, the same value that the zero vector gives.
curve <- quadratic_penalty(diag(2), map = diff(diag(4), differences = 2))
penalty_value(curve, c(1, 2, 3, 4), list(lambda = 5))
#> [1] -1.269427
penalty_value(curve, c(0, 0, 0, 0), list(lambda = 5))
#> [1] -1.269427
penalty_value(curve, c(1, 2, 4, 8), list(lambda = 5))
#> [1] 11.23057
```
