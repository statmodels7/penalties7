# Construct a Structured Quadratic Penalty

Builds the exact negative log-density of the Gaussian prior \\\beta \sim
N(0, \Omega(\theta)^{-1})\\, where the precision \\\Omega(\theta)\\ is a
parameters7 matrix parameter and the hyperparameters are its free
values. This is the correlated prior that
[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
cannot express: in a quadratic penalty one scale multiplies a constant
matrix, while here the hyperparameters reach every entry of the matrix.

## Usage

``` r
structured_penalty(structure)
```

## Arguments

- structure:

  A parameters7 `matrix_parameter`, read as the prior's precision. Pass
  `parameters7::inverse_of(s)` to put `s` on the covariance side
  instead.

  It may be rank deficient, giving an improper prior:
  [`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
  then returns `FALSE` and the constant uses the rank and the log
  pseudo-determinant. A rank-deficient structure can be used only on the
  precision side, because
  [`parameters7::inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.html)
  rejects it, a singular matrix having no inverse.

## Value

A
[`StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/StructuredPenalty.md)
object whose hyperparameters are the structure's `free_names`, each
unconstrained and on an identity link.

## The value

Writing \\\Omega\\ for the precision, which is the matrix that the
structure supplies,

\$\$\rho(\beta; \theta) = \tfrac{1}{2}\\\beta'\Omega(\theta)\beta -
\tfrac{1}{2}\log\mathrm{pdet}\\\Omega(\theta) + \tfrac{r}{2}\log
2\pi,\$\$

with \\r\\ the structure's rank. For a full-rank structure this is the
negative log-density of the multivariate normal distribution with mean
zero and precision \\\Omega(\theta)\\; the examples compare it with a
direct evaluation of that density.

## The structure is the precision

The matrix that the structure supplies is the prior's **precision**:
\\\Omega = M(\eta)\\. This follows from the construction and is not an
option, because \\\rho\\ is a negative log-density and the matrix that
enters its quadratic form is the precision.

The choice of the side on which a structure is placed still matters,
because the two readings are different priors. For a positive definite
matrix they agree only where \\M = M^{-1}\\, that is at the identity.
The inverse of an AR(1) matrix is tridiagonal, so for a dimension of
three or more it is an AR(1) matrix only at zero correlation. To place a
structure `s` on the covariance side, build the family whose value is
its inverse,

    structured_penalty(parameters7::inverse_of(s))

which carries \\\partial_k \Sigma^{-1} = -\Sigma^{-1}B_k\Sigma^{-1}\\,
with \\B_k = \partial_k\Sigma\\ the derivative of the covariance, and
its higher orders inside the structure, where the log-determinant stays
exact. For
[`parameters7::ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.html)
and
[`parameters7::autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.html)
the inverse has a constructor of its own,
[`parameters7::ar1_inv()`](https://statmodels7.github.io/parameters7/reference/ar1_inv.html)
and
[`parameters7::autoregressive_inv()`](https://statmodels7.github.io/parameters7/reference/autoregressive_inv.html),
and a structure closed under inversion requires no wrapper.

## The derivatives

Every derivative is computed from the derivative arrays that the
structure supplies. The hyperparameter gradient is
\\\tfrac{1}{2}\beta'A_k\beta -
\tfrac{1}{2}\partial_k\log\mathrm{pdet}\\\Omega\\, with \\A_k =
\partial_k\Omega\\ the structure's `param_d1`; the Hessian uses
`param_d2`, and the mixed block is \\A_k\beta\\. No conversion between
covariance and precision happens in this package: a structure meant for
the covariance side is a different structure, and it supplies its own
arrays.

## No map

There is no `map` argument. A linear image of a structured precision is
a different precision, and composing it into the structure, where its
log-determinant stays exact, is left to parameters7.

## See also

[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
for one scale on a fixed matrix,
[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
for a sum of quadratics,
[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md)
for a coordinatewise prior,
[`parameters7::log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.html)
and
[`parameters7::ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.html)
for structures to pass in,
[`penalty_readable()`](https://statmodels7.github.io/penalties7/reference/penalty_readable.md)
for the interpretable quantities (a scale, a correlation) of a structure
that declares them through
[`parameters7::param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.html).

## Examples

``` r
# An AR(1) prior on four coefficients: two hyperparameters reach every
# entry of the precision.
pen <- structured_penalty(parameters7::ar1(4))
pen@params
#> [1] "log_scale" "z_rho"    
theta <- list(log_scale = 0.2, z_rho = 0.5)
b <- c(0.3, -0.1, 0.4, 0.2)
penalty_value(pen, b, theta)
#> [1] 3.858268

# The Hessian is the matrix of the structure, read as a precision.
Om <- parameters7::param_value(parameters7::ar1(4),
                               c(0.2, 0.5))
max(abs(penalty_hessian(pen, b, theta) - unclass(Om)))
#> [1] 0

# The same structure read as a covariance is a different prior at the same
# free values: its Hessian is the inverse of that matrix.
cov <- structured_penalty(parameters7::inverse_of(parameters7::ar1(4)))
max(abs(penalty_hessian(cov, b, theta) - solve(unclass(Om))))
#> [1] 5.551115e-17
c(precision = penalty_value(pen, b, theta),
  covariance = penalty_value(cov, b, theta))
#>  precision covariance 
#>   3.858268   3.885654 

# The covariance reading is the negative log-density of the multivariate
# normal with covariance S, written out directly.
S <- unclass(Om)
penalty_value(cov, b, theta) -
  (0.5 * sum(b * solve(S, b)) +
     0.5 * as.numeric(determinant(S, logarithm = TRUE)$modulus) +
     2 * log(2 * pi))
#> [1] 0

# At a zero log-Cholesky free vector the structure is the identity, so the
# penalty equals the plain ridge at lambda = 1, up to rounding.
s <- structured_penalty(parameters7::log_cholesky(3))
z <- as.list(stats::setNames(rep(0, length(s@params)), s@params))
bb <- c(0.4, -1.1, 0.7)
penalty_value(s, bb, z) - penalty_value(ridge_penalty(n_coef = 3), bb,
                                        list(lambda = 1))
#> [1] 0
```
