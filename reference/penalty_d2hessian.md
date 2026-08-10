# The Second Derivative of the Coefficient Hessian in the Hyperparameters

\\\partial^2 S/\partial\theta_m\partial\theta_l =
\partial^4\rho/\partial\beta^2\partial\theta_m\partial\theta_l\\, one
matrix per unordered pair.

## Usage

``` r
penalty_d2hessian(pen, beta, theta, ...)
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

A named list of matrices, keyed by hyperparameter pair.

## Details

Zero for every penalty whose Hessian is linear in its hyperparameters,
which is the quadratic and additive branches; the structured branch
reads the matrix parameter's `param_d2`, and the separable branch the
second \\\theta\\-derivative of the parent's response curvature.

The keys are those of
[`penalty_hess_theta`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md),
so a consumer looking a pair up need not know which order it was written
in.

## See also

[`penalty_dhessian`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)

## Examples

``` r
pen <- quadratic_penalty(diag(3))
penalty_d2hessian(pen, c(1, 2, 3), list(lambda = 2))
#> $lambda_lambda
#>      [,1] [,2] [,3]
#> [1,]    0    0    0
#> [2,]    0    0    0
#> [3,]    0    0    0
#> 
```
