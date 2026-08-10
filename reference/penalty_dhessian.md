# The Derivative of the Coefficient Hessian in the Hyperparameters

\\\partial S/\partial\theta_m = \partial^3\rho/\partial\beta^2\\
\partial\theta_m\\, one matrix per hyperparameter.

## Usage

``` r
penalty_dhessian(pen, beta, theta, ...)
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

A named list of matrices, one per hyperparameter.

## Details

`penalty_hessian` says how curved the penalty is in the coefficients;
this says how that curvature moves with each hyperparameter. It is the
piece a marginal criterion needs to differentiate \\\log\|H+S\|\\.

Every branch answers in closed form. For
[`quadratic_penalty`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
the Hessian is \\\lambda D'PD\\ and the derivative is \\D'PD\\; for
[`additive_penalty`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
it is the component \\P_k\\ of each smoothing parameter; for
[`structured_penalty`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md)
it is the matrix parameter's own `param_d1`; for
[`distrib_penalty`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md)
it is \\-D'\mathrm{diag}(\partial^3\ell/\partial
y^2\partial\theta_m)D\\, which distributions7 supplies as
`distrib_cross2_y`. A penalty defined by its derivative rather than by a
density –
[`scad_penalty`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md),
[`mcp_penalty`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
– has no such quantity where its kinks are and rejects.

## See also

[`penalty_hessian`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md),
[`penalty_d2hessian`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md),
[`penalty_dcross`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md)

## Examples

``` r
pen <- quadratic_penalty(diag(3))
penalty_dhessian(pen, c(1, 2, 3), list(lambda = 2))
#> $lambda
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
#> 
```
