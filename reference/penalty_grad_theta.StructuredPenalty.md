# Hyperparameter Derivatives of a Structured Penalty

The three blocks in the structure's free values, all assembled from the
structure's own derivative arrays.
[`penalty_grad_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
returns one number per free value,
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
one per unordered pair, and
[`penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
one coefficient vector per free value.

## Arguments

- pen:

  A
  [`StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/StructuredPenalty.md)
  object.

- beta:

  A numeric vector of length `pen@n_coef`.

- theta:

  A named list of the structure's free values.

- scale:

  Read by the generic, which applies the chain rule after dispatch.
  Every link here is the identity, so the two scales coincide and
  `"link"` returns the same numbers.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`penalty_grad_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
a list of one number per free value, named by `pen@params`.
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
a list of one number per unordered pair, keyed diagonals first.
[`penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
a list of one numeric vector of length `pen@n_coef` per free value,
named by `pen@params`.

## Details

With \\A_k = \partial_k\Omega\\ and \\A\_{kl} = \partial\_{kl}\Omega\\
the precision's derivative arrays, transported from the covariance where
the structure describes one,

\$\$\frac{\partial\rho}{\partial\theta_k} = \tfrac{1}{2}\beta'A_k\beta -
\tfrac{1}{2}\partial_k \log\mathrm{pdet}\\\Omega, \qquad
\frac{\partial^2\rho}{\partial\theta_k\partial\theta_l} =
\tfrac{1}{2}\beta'A\_{kl}\beta - \tfrac{1}{2}\partial\_{kl}
\log\mathrm{pdet}\\\Omega, \qquad
\frac{\partial^2\rho}{\partial\beta\\\partial\theta_k} = A_k\beta.\$\$

The mixed block carries no log-determinant term, that term not depending
on the coefficients.

The second-order components are keyed by pairs as
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
keys them, diagonals first and then the upper off-diagonal pairs joined
by an underscore. Internally they are looked up in the structure's own
keying, the two free names joined by a colon and sorted by position, and
the key is constructed from the pair rather than parsed out of a name.

## See also

[`penalty_value.StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.StructuredPenalty.md)
for the quantity differentiated,
[`struct_omega()`](https://statmodels7.github.io/penalties7/reference/struct_omega.md)
for the arrays these read,
[`penalty_logpdet.StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.StructuredPenalty.md)
for the determinant's own derivatives.

## Examples

``` r
pen <- structured_penalty(parameters7::ar1(4, role = "precision"))
th <- list(log_scale = 0.2, z_rho = 0.5)
b <- c(0.3, -0.1, 0.4, 0.2)

penalty_grad_theta(pen, b, th)
#> $log_scale
#> [1] -1.77783
#> 
#> $z_rho
#> [1] 1.52166
#> 
names(penalty_hess_theta(pen, b, th))
#> [1] "log_scale_log_scale" "z_rho_z_rho"         "log_scale_z_rho"    

# Every link is the identity, so the unconstrained scale is the same one.
all.equal(penalty_grad_theta(pen, b, th),
          penalty_grad_theta(pen, b, th, scale = "link"))
#> [1] TRUE

# The mixed block is A_k beta, so summing it against beta gives twice the
# quadratic part of the gradient in that direction.
cr <- penalty_cross(pen, b, th)
sum(b * cr$z_rho) / 2
#> [1] 0.1353085
```
