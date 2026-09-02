# Construct a Structured Quadratic Penalty

Builds the exact negative log-density of the Gaussian prior \\\beta \sim
N(0, \Sigma(\theta))\\, where the matrix is a parameters7 matrix
parameter and the hyperparameters are its free values. This is the
correlated prior
[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
cannot express: there one scale multiplies a constant matrix, here the
hyperparameters reach every entry.

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
  then answers `FALSE` and the constant uses the rank and the log
  pseudo-determinant. That is the one case the covariance reading cannot
  have, and
  [`parameters7::inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.html)
  rejects it there for the reason a deficient covariance has no inverse.

## Value

A
[`StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/StructuredPenalty.md)
object whose hyperparameters are the structure's `free_names`, each
unconstrained and on an identity link.

## The value

Writing \\\Omega\\ for the precision, which is what the structure
supplies,

\$\$\rho(\beta; \theta) = \tfrac{1}{2}\\\beta'\Omega(\theta)\beta -
\tfrac{1}{2}\log\mathrm{pdet}\\\Omega(\theta) + \tfrac{r}{2}\log
2\pi,\$\$

with \\r\\ the structure's rank. For a full-rank structure this is
exactly `-mvtnorm`-style multivariate normal log-density, which the
examples check against a hand-written one.

## The structure is the precision

The matrix the structure supplies is the prior's **precision**: \\\Omega
= M(\eta)\\. That is a property of the construction and not a choice
offered to the caller, \\\rho\\ being a negative log-density and the
matrix it needs being the one in the quadratic form.

Which side a structure is meant for still matters, because the two are
different priors: there is in general no free vector at which they
agree, the inverse of an AR(1) covariance being tridiagonal and not an
AR(1) at any parameters. A caller who wants the structure on the
covariance side says so by building the family whose value is its
inverse,

    structured_penalty(parameters7::inverse_of(s))

which carries \\\partial_k \Sigma^{-1} = -\Sigma^{-1}A_k\Sigma^{-1}\\
and its higher orders inside the structure, where the log-determinant
stays exact. For
[`parameters7::ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.html)
and
[`parameters7::autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.html)
the inverse has a name of its own,
[`parameters7::ar1_inv()`](https://statmodels7.github.io/parameters7/reference/ar1_inv.html)
and
[`parameters7::autoregressive_inv()`](https://statmodels7.github.io/parameters7/reference/autoregressive_inv.html),
and for a structure closed under inversion nothing is needed at all.

## The derivatives

Every derivative comes from the structure's own contract. The
hyperparameter gradient is \\\tfrac{1}{2}\beta'A_k\beta -
\tfrac{1}{2}\partial_k\log\mathrm{pdet}\\ with \\A_k\\ the structure's
`param_d1`, the Hessian adds `param_d2`, and the mixed block is
\\A_k\beta\\. Nothing is transported here: a structure meant for the
other side is a different structure, and it supplies its own arrays.

## No map

There is no `map` argument. A linear image of a structured precision is
a different precision, and composing it into the structure, where its
log-determinant stays exact, is the structure's own business.

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
for reporting the hyperparameters as standard deviations and
correlations.

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

# The Hessian is the structure's own matrix, read as a precision.
Om <- parameters7::param_value(parameters7::ar1(4),
                               c(0.2, 0.5))
max(abs(penalty_hessian(pen, b, theta) - unclass(Om)))
#> [1] 0

# The same structure read as a covariance is a DIFFERENT prior at the same
# free values: there the Hessian is the inverse of that matrix.
cov <- structured_penalty(parameters7::inverse_of(parameters7::ar1(4)))
max(abs(penalty_hessian(cov, b, theta) - solve(unclass(Om))))
#> [1] 1.110223e-16
c(precision = penalty_value(pen, b, theta),
  covariance = penalty_value(cov, b, theta))
#>  precision covariance 
#>   3.858268   3.885654 

# And the covariance reading is exactly a multivariate normal log-density.
S <- unclass(Om)
penalty_value(cov, b, theta) -
  (0.5 * sum(b * solve(S, b)) +
     0.5 * as.numeric(determinant(S, logarithm = TRUE)$modulus) +
     2 * log(2 * pi))
#> [1] 0

# At a zero log-Cholesky free vector the structure is the identity, so the
# penalty is the plain ridge at lambda = 1, to the last bit.
s <- structured_penalty(parameters7::log_cholesky(3))
z <- as.list(stats::setNames(rep(0, length(s@params)), s@params))
bb <- c(0.4, -1.1, 0.7)
penalty_value(s, bb, z) - penalty_value(ridge_penalty(n_coef = 3), bb,
                                        list(lambda = 1))
#> [1] 0
```
