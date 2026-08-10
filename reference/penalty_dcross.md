# The Derivative of the Mixed Block in the Hyperparameters

\\\partial^3\rho/\partial\beta\\\partial\theta_m\partial\theta_l\\, one
vector per unordered pair.

## Usage

``` r
penalty_dcross(pen, beta, theta, ...)
```

## Arguments

- pen:

  A
  [`penalty`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- beta:

  A numeric vector of coefficients.

- theta:

  A named list of hyperparameters.

- ...:

  Passed to methods.

## Value

A named list of numeric vectors, keyed by hyperparameter pair.

## Details

[`penalty_cross`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
is how the coefficient gradient moves with one hyperparameter; this is
how that movement itself moves with a second. It is what a marginal
criterion needs to differentiate the mode's own derivative, and it is
exactly zero wherever the penalty is quadratic in the coefficients with
a Hessian linear in the hyperparameters.

## See also

[`penalty_cross`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md),
[`penalty_dhessian`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)

## Examples

``` r
pen <- quadratic_penalty(diag(3))
penalty_dcross(pen, c(1, 2, 3), list(lambda = 2))
#> $lambda_lambda
#> [1] 0 0 0
#> 
```
