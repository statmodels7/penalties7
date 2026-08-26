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

  A parameters7 `matrix_parameter` whose `role` says which matrix of the
  prior it is.

  A structure declared `"precision"` may be rank deficient, giving an
  improper prior:
  [`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
  then answers `FALSE` and the constant uses the rank and the log
  pseudo-determinant.

  A structure declared `"covariance"` may not be rank deficient. A
  covariance of deficient rank has no inverse, and a direction of zero
  variance is a constraint on the coefficients and not a prior over
  them; the constructor rejects it, naming the rank and the dimension.

  A structure declaring `"either"` is rejected, with the two roles
  named.

## Value

A
[`StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/StructuredPenalty.md)
object whose hyperparameters are the structure's `free_names`, each
unconstrained and on an identity link.

## The value

Writing \\\Omega\\ for the precision, whether the structure supplies it
directly or through the covariance it describes,

\$\$\rho(\beta; \theta) = \tfrac{1}{2}\\\beta'\Omega(\theta)\beta -
\tfrac{1}{2}\log\mathrm{pdet}\\\Omega(\theta) + \tfrac{r}{2}\log
2\pi,\$\$

with \\r\\ the structure's rank. For a full-rank structure this is
exactly `-mvtnorm`-style multivariate normal log-density, which the
examples check against a hand-written one.

## The two roles

The structure's `role` says which matrix of the prior it is, and the two
readings are different priors from the same free vector. A structure of
role `"precision"` at \\\eta\\ gives \\\Omega = M(\eta)\\; one of role
`"covariance"` gives \\\Omega = M(\eta)^{-1}\\. There is in general no
free vector that makes the two agree: the inverse of an AR(1) covariance
is tridiagonal and is not an AR(1) covariance at any parameters.

A structure that declares `"either"` is rejected. The two readings
differ in the sign of the log-determinant term, and nothing in the
matrix says which was meant, so a default would give a fit that
converges to a different prior without saying so.

## The derivatives

Every derivative comes from the structure's own contract. Where the
structure is the precision, the hyperparameter gradient is
\\\tfrac{1}{2}\beta'A_k\beta - \tfrac{1}{2}\partial_k\log\mathrm{pdet}\\
with \\A_k\\ the structure's `param_d1`, the Hessian adds `param_d2`,
and the mixed block is \\A_k\beta\\. Where it is the covariance the same
expressions are read at the precision it implies, whose derivatives
follow from the chain rule for an inverse,

\$\$\partial_k\Omega = -\Omega A_k \Omega, \qquad \partial\_{kl}\Omega =
\Omega\left(A_k\Omega A_l + A_l\Omega A_k\right)\Omega - \Omega
A\_{kl}\Omega,\$\$

with \\\log\lvert\Omega\rvert = -\log\lvert\Sigma\rvert\\ and its
derivatives negated termwise. The transport is done once, in
[`struct_omega()`](https://statmodels7.github.io/penalties7/reference/struct_omega.md)
and its siblings, so the methods are the same arithmetic in both cases.

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
pen <- structured_penalty(parameters7::ar1(4, role = "precision"))
pen@params
#> [1] "log_scale" "z_rho"    
theta <- list(log_scale = 0.2, z_rho = 0.5)
b <- c(0.3, -0.1, 0.4, 0.2)
penalty_value(pen, b, theta)
#> [1] 3.858268

# The Hessian is the structure's own matrix, read as a precision.
Om <- parameters7::param_value(parameters7::ar1(4, role = "precision"),
                               c(0.2, 0.5))
max(abs(penalty_hessian(pen, b, theta) - unclass(Om)))
#> [1] 0

# The same structure read as a covariance is a DIFFERENT prior at the same
# free values: there the Hessian is the inverse of that matrix.
cov <- structured_penalty(parameters7::ar1(4, role = "covariance"))
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
s <- structured_penalty(parameters7::log_cholesky(3, role = "precision"))
z <- as.list(stats::setNames(rep(0, length(s@params)), s@params))
bb <- c(0.4, -1.1, 0.7)
penalty_value(s, bb, z) - penalty_value(ridge_penalty(n_coef = 3), bb,
                                        list(lambda = 1))
#> [1] 0
```
