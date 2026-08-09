# The Non-Differentiable Points of a Penalty

The values of \\t = D\beta\\ at which \\\rho\\ is not differentiable in
its argument: empty for the smooth penalties, `0` for the lasso, SCAD
and MCP.
[`check_penalty`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
places its grids away from them, and a solver may consult them.

## Usage

``` r
penalty_kinks(pen, theta, ...)
```

## Arguments

- pen:

  A
  [`penalty`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- theta:

  A named list of hyperparameter values.

- ...:

  Passed to methods.

## Value

A numeric vector, possibly empty.

## See also

[`penalty_value`](https://statmodels7.github.io/penalties7/reference/penalty_value.md),
[`penalty_gradient`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md),
[`penalty_hessian`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md),
[`penalty_grad_theta`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md),
[`penalty_cross`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)

## Examples

``` r
penalty_kinks(quadratic_penalty(diag(2)), list(lambda = 1))
#> numeric(0)
penalty_kinks(lasso_penalty(), list(lambda = 1))
#> [1] 0
```
