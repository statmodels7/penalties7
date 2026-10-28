# The Derivative of the Mixed Block in the Hyperparameters

Returns
\\\partial^3\rho/\partial\beta\\\partial\theta_m\partial\theta_l\\, one
coefficient vector per unordered pair.
[`penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
gives the derivative of the coefficient gradient in one hyperparameter,
and this generic gives the derivative of that quantity in a second
hyperparameter.

## Usage

``` r
penalty_dcross(penalty, beta, theta, ...)
```

## Arguments

- penalty:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- beta:

  A numeric vector of length `penalty@n_coef`; the generic signals an
  error for any other length. The structured and separable branches read
  its values; the quadratic and additive branches return zero for any
  value.

- theta:

  A named list of hyperparameter values, or a named numeric vector
  carrying the same, checked by the generic as in
  [`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md).

- ...:

  Passed to methods. The shipped methods do not read it.

## Value

A named list of one numeric vector of length `penalty@n_coef` per
unordered pair, keyed diagonals first.

## Details

A marginal criterion is evaluated at the penalized mode, and the mode
moves with the hyperparameters. The second derivative of the criterion
therefore contains the second derivative of the mode, and this quantity
is the contribution of the penalty to it.

It is exactly zero wherever the penalty is quadratic in the coefficients
with a Hessian linear in the hyperparameters, which covers the quadratic
and additive branches. The structured branch returns
\\(\partial^2\Omega/\partial\theta_m\partial\theta_l)\beta\\, computed
from
[`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md),
and the separable branch returns \\-D'(\partial^3\ell/\partial
y\\\partial\theta_m\partial\theta_l)\\, from
[`distributions7::distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.html)
of the parent evaluated at \\D\beta\\.

The keys are those of
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md),
diagonal pairs first.

## See also

[`penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
for the first order,
[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
and
[`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md)
for the other two quantities that a marginal criterion requires.

## Examples

``` r
b <- c(1, -0.5, 0.3)

# The derivative is zero on the quadratic branch.
penalty_dcross(quadratic_penalty(diag(3)), b, list(lambda = 2))
#> $lambda_lambda
#> [1] 0 0 0
#> 

# It is not zero on the structured branch, where it equals
# d2 Omega / dtheta^2 times the coefficients.
s <- structured_penalty(parameters7::log_cholesky(2))
th <- list(log_L1 = 0.1, log_L2 = -0.1, L2.1 = 0.3)
bb <- c(1, -0.5)
penalty_dcross(s, bb, th)$log_L1_log_L1
#> [1] 4.7198354 0.3315513
drop(penalty_d2hessian(s, bb, th)$log_L1_log_L1 %*% bb)
#> [1] 4.7198354 0.3315513
```
