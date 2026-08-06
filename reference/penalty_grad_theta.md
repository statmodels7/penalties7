# Hyperparameter Derivatives of a Penalty

`penalty_grad_theta` returns \\\partial\rho/\partial\theta\\ as a named
list, `penalty_hess_theta` the second derivatives keyed with diagonals
first, and `penalty_cross` the mixed block
\\\partial^2\rho/\partial\beta\\\partial\theta_k\\, one coefficient
vector per hyperparameter – the block a joint estimation of coefficients
and hyperparameters needs.

## Usage

``` r
penalty_grad_theta(pen, beta, theta, scale = c("parameter", "link"), ...)

penalty_hess_theta(pen, beta, theta, scale = c("parameter", "link"), ...)

penalty_cross(pen, beta, theta, scale = c("parameter", "link"), ...)
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

- scale:

  Either `"parameter"` (default) or `"link"`; on the link scale the
  derivatives are with respect to the unconstrained values, carried by
  the chain rule in the generic body, so methods always return the
  parameter scale.

- ...:

  Passed to methods.

## Value

A named list: one number per hyperparameter for the gradient, one number
per pair for the Hessian, one numeric vector of length `q` per
hyperparameter for the mixed block.

## Examples

``` r
pen <- quadratic_penalty(diag(2))
penalty_grad_theta(pen, c(1, -1), list(lambda = 3))
#> $lambda
#> [1] 0.6666667
#> 
penalty_cross(pen, c(1, -1), list(lambda = 3))
#> $lambda
#> [1]  1 -1
#> 
```
