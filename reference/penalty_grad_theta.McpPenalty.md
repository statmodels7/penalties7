# Hyperparameter Derivatives of an MCP Penalty

The three blocks in \\(\lambda, \gamma)\\, each a closed piecewise form
summed over the coordinates.
[`penalty_grad_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
returns the two first derivatives,
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
the three second derivatives, and
[`penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
the mixed block, one coefficient vector per hyperparameter.

## Arguments

- pen:

  An
  [`McpPenalty()`](https://statmodels7.github.io/penalties7/reference/ScadPenalty.md)
  object.

- beta:

  A numeric vector of length `pen@n_coef`.

- theta:

  A list holding `lambda` and `gamma`.

- scale:

  Read by the generic, which applies the chain rule onto the
  unconstrained scale after dispatch. These methods always return the
  parameter scale.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`penalty_grad_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
a list of two numbers, `lambda` and `gamma`.
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
a list of three numbers, `lambda_lambda`, `gamma_gamma` and
`lambda_gamma`.
[`penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
a list of two numeric vectors of length `pen@n_coef`, named `lambda` and
`gamma`.

## Details

Inside the taper the value is \\\lambda u - u^2/(2\gamma)\\, so
\\\partial\rho/\partial\lambda = u\\ and \\\partial\rho/\partial\gamma =
u^2/(2\gamma^2)\\; beyond it the value is \\\gamma\lambda^2/2\\ and the
two derivatives are \\\gamma\lambda\\ and \\\lambda^2/2\\, neither
depending on the coefficient. The second derivative in \\\lambda\\ is
therefore zero inside the taper, the value being linear in \\\lambda\\
there.

As for SCAD, these derivatives exist and no marginal criterion can use
them: the expansion such a criterion rests on needs a mode at which the
penalty is twice differentiable, and a mode with a coefficient at zero
is on the kink.

## See also

[`penalty_value.McpPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.McpPenalty.md)
for the quantity differentiated,
[`penalty_grad_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
for the generic and its two scales,
[`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md),
which verifies all three.

## Examples

``` r
pen <- mcp_penalty(n_coef = 3)
th <- list(lambda = 1, gamma = 3)
b <- c(0.5, 2, 5)

penalty_grad_theta(pen, b, th)
#> $lambda
#> [1] 5.5
#> 
#> $gamma
#> [1] 0.7361111
#> 

# The value is linear in lambda inside the taper, so only the coordinate
# beyond gamma * lambda contributes to the second derivative.
penalty_hess_theta(pen, b, th)$lambda_lambda
#> [1] 3
```
