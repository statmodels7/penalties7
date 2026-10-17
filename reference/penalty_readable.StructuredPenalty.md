# Interpretable Quantities of a Structured Penalty

Returns the interpretable quantities that the structure declares through
[`parameters7::param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.html),
in the shape that
[`penalty_readable.DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_readable.DistribPenalty.md)
returns, and `NULL` for a structure that declares none.

## Arguments

- pen:

  A
  [`StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/StructuredPenalty.md)
  object.

- theta:

  A named list of the structure's free values.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

`NULL` when the structure declares no quantities. Otherwise the function
returns a list of four: `value`, a numeric vector of the quantities, in
which each name is the structure's own name with the prefix
`precision_`, or `covariance_` for an inverse structure; `jacobian`,
their derivatives, one row per quantity and one column per
hyperparameter, with row names from `value` and column names from
`pen@params`; `transform`, a character vector named like `value` that
names the scale on which each interval is built; and `block`, the
structure's label for every quantity, named like `value`.

## Details

The hyperparameters of a structured penalty are the free values of its
structure, such as the logarithm of a scale and the inverse hyperbolic
tangent of a correlation for
[`parameters7::ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.html).
[`parameters7::param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.html)
maps them to the quantities that the structure declares. Among the
shipped structures,
[`parameters7::ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.html)
and
[`parameters7::compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.html)
declare a scale and a correlation,
[`parameters7::scaled_matrix()`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.html)
its multiplier, and
[`parameters7::autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.html)
the scale, the partial autocorrelations and the autoregressive
coefficients.
[`parameters7::inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.html),
[`parameters7::ar1_inv()`](https://statmodels7.github.io/parameters7/reference/ar1_inv.html)
and
[`parameters7::autoregressive_inv()`](https://statmodels7.github.io/parameters7/reference/autoregressive_inv.html)
declare the quantities, if any, of the structure that they invert. The
other structures, among them
[`parameters7::log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.html)
and the composition wrappers, declare nothing, and the method then
returns `NULL`.

The quantities describe the matrix of the structure, and for a structure
that is not an inverse that matrix is the prior's precision. For
`structured_penalty(parameters7::ar1(4))` the scale is the diagonal
entry of the precision and the correlation is the ratio of an entry next
to the diagonal to the diagonal entry; they are not the variance and the
correlation of the coefficients. For this reason the method prefixes
every name that
[`parameters7::param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.html)
gives with `precision_`, so that `scale` and `rho` are reported as
`precision_scale` and `precision_rho`.

For an inverse structure
[`parameters7::param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.html)
marks the declaration with `inverse = TRUE`, and the quantities describe
the inverse of the precision, which is the prior's covariance. For
`structured_penalty(parameters7::inverse_of(parameters7::ar1(4)))` they
are the variance and the lag-one correlation of the coefficients, and
the method prefixes them with `covariance_`, so that they are reported
as `covariance_scale` and `covariance_rho`.

Every link of the branch is the identity, so the hyperparameters are the
free values themselves, and the Jacobian that
[`parameters7::param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.html)
computes in closed form is the Jacobian with respect to the
hyperparameters. The method names its rows by the prefixed names and its
columns by `pen@params`, and repeats the structure's label once per
quantity to form `block`.

## See also

[`penalty_readable()`](https://statmodels7.github.io/penalties7/reference/penalty_readable.md)
for the generic,
[`parameters7::param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.html)
for the quantities that a structure declares.

## Examples

``` r
# The scale and the correlation of an AR(1) precision on four coefficients
# are named precision_scale and precision_rho.
pen <- structured_penalty(parameters7::ar1(4))
r <- penalty_readable(pen, list(log_scale = 0.2, z_rho = 0.5))
r$value
#> precision_scale   precision_rho 
#>       1.2214028       0.4621172 
r$jacobian
#>                 log_scale     z_rho
#> precision_scale  1.221403 0.0000000
#> precision_rho    0.000000 0.7864477

# The quantities describe the precision: its diagonal entry, and the ratio
# of an entry next to the diagonal to the diagonal entry.
Om <- penalty_matrix(pen, list(log_scale = 0.2, z_rho = 0.5))
c(Om[1, 1], Om[1, 2] / Om[1, 1])
#> [1] 1.2214028 0.4621172

# On the covariance side the same free values give the variance and the
# correlation of the coefficients, named covariance_scale and
# covariance_rho.
cov <- structured_penalty(parameters7::inverse_of(parameters7::ar1(4)))
rc <- penalty_readable(cov, list(log_scale = 0.2, z_rho = 0.5))
rc$value
#> covariance_scale   covariance_rho 
#>        1.2214028        0.4621172 

# They are read off the inverse of the precision.
Sigma <- solve(penalty_matrix(cov, list(log_scale = 0.2, z_rho = 0.5)))
c(Sigma[1, 1], Sigma[1, 2] / Sigma[1, 1])
#> [1] 1.2214028 0.4621172

# A log-Cholesky structure declares nothing.
lc <- structured_penalty(parameters7::log_cholesky(2))
penalty_readable(lc, list(log_L1 = 0, log_L2 = 0, L2.1 = 0.3))
#> NULL
```
