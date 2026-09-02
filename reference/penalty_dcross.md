# The Derivative of the Mixed Block in the Hyperparameters

Returns
\\\partial^3\rho/\partial\beta\\\partial\theta_m\partial\theta_l\\, one
coefficient vector per unordered pair. Where
[`penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
says how the coefficient gradient moves with one hyperparameter, this
says how that movement itself moves with a second.

## Usage

``` r
penalty_dcross(pen, beta, theta, ...)
```

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- beta:

  A numeric vector of length `pen@n_coef`. Read by the structured and
  separable branches; the quadratic and additive ones answer zero
  whatever it is.

- theta:

  A named list of hyperparameter values, or a named numeric vector
  carrying the same.

- ...:

  Passed to methods. No shipped method reads it.

## Value

A named list of one numeric vector of length `pen@n_coef` per unordered
pair, keyed diagonals first.

## Details

A marginal criterion is evaluated at the penalized mode, and the mode
moves with the hyperparameters. Differentiating the criterion twice
therefore needs the derivative of the mode's own derivative, and this is
the penalty's contribution to it.

It is **exactly zero** wherever the penalty is quadratic in the
coefficients with a Hessian linear in the hyperparameters, which covers
the quadratic and additive branches. The structured branch answers
\\(\partial^2\Omega/\partial\theta_m\partial\theta_l)\beta\\, reusing
[`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md),
and the separable branch reads the parent's
[`distributions7::distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.html).

The keys are
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)'s,
diagonals first.

## See also

[`penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
for the first order,
[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
and
[`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md)
for the other two quantities a marginal criterion asks for.

## Examples

``` r
b <- c(1, -0.5, 0.3)

# Zero on the quadratic branch.
penalty_dcross(quadratic_penalty(diag(3)), b, list(lambda = 2))
#> $lambda_lambda
#> [1] 0 0 0
#> 

# And not on the structured one, where it is d2 Omega / dtheta^2 times
# the coefficients.
s <- structured_penalty(parameters7::log_cholesky(2))
th <- list(log_L1 = 0.1, log_L2 = -0.1, L2.1 = 0.3)
bb <- c(1, -0.5)
penalty_dcross(s, bb, th)$log_L1_log_L1
#> [1] 4.7198354 0.3315513
drop(penalty_d2hessian(s, bb, th)$log_L1_log_L1 %*% bb)
#> [1] 4.7198354 0.3315513
```
