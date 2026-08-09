# Coefficient Derivatives of a Penalty

`penalty_gradient` returns \\\partial\rho/\partial\beta\\ and
`penalty_hessian` returns \\\partial^2\rho/\partial\beta^2\\, both
exact.

## Usage

``` r
penalty_gradient(pen, beta, theta, ...)

penalty_hessian(pen, beta, theta, ...)
```

## Arguments

- pen:

  A
  [`penalty`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- beta:

  A numeric vector of coefficients.

- theta:

  A named list of hyperparameter values.

- ...:

  Passed to methods.

## Value

`penalty_gradient` a numeric vector of length `q`; `penalty_hessian` a
`q x q` symmetric matrix.

## See also

[`penalty_value`](https://statmodels7.github.io/penalties7/reference/penalty_value.md),
[`penalty_grad_theta`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md),
[`penalty_kinks`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)

## Examples

``` r
pen <- quadratic_penalty(diag(2))
penalty_gradient(pen, c(1, -1), list(lambda = 3))
#> [1]  3 -3
penalty_hessian(pen, c(1, -1), list(lambda = 3))
#>      [,1] [,2]
#> [1,]    3    0
#> [2,]    0    3
```
