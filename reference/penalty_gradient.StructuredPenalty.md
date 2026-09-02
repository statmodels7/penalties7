# Coefficient Derivatives of a Structured Penalty

[`penalty_gradient()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
returns \\\Omega(\theta)\beta\\ and
[`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
returns \\\Omega(\theta)\\ itself. The value is a quadratic form in
\\\beta\\, so the Hessian does not move with the coefficients and the
gradient is the Hessian applied to them.

## Arguments

- pen:

  A
  [`StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/StructuredPenalty.md)
  object.

- beta:

  A numeric vector of length `pen@n_coef`.

- theta:

  A named list of the structure's free values.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`penalty_gradient()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
an unnamed numeric vector of length `pen@n_coef`;
[`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
the symmetric precision matrix of that side, **carrying the structure's
own dimnames** (`v1`, `v2`, ... for the parameters7 primitives). The
other branches return matrices with no dimnames, so a comparison across
branches wants [`unname()`](https://rdrr.io/r/base/unname.html).

## Details

Where the structure is a covariance, both read the inverse it implies,
so each call costs one solve. A caller taking several derivatives at one
\\\theta\\ pays for that inversion once per generic.

## See also

[`penalty_value.StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.StructuredPenalty.md)
for the quantity differentiated,
[`penalty_grad_theta.StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.StructuredPenalty.md)
for the hyperparameter blocks.

## Examples

``` r
pen <- structured_penalty(parameters7::ar1(4))
th <- list(log_scale = 0.2, z_rho = 0.5)
b <- c(0.3, -0.1, 0.4, 0.2)

# The gradient is the constant Hessian applied to the coefficients. The
# matrix carries the structure's dimnames and the vector does not.
all.equal(penalty_gradient(pen, b, th),
          unname(drop(penalty_hessian(pen, b, th) %*% b)))
#> [1] TRUE
dimnames(penalty_hessian(pen, b, th))[[1]]
#> [1] "v1" "v2" "v3" "v4"

# And the Hessian is the structure's own matrix.
Om <- parameters7::param_value(parameters7::ar1(4),
                               c(0.2, 0.5))
max(abs(penalty_hessian(pen, b, th) - unclass(Om)))
#> [1] 0
```
