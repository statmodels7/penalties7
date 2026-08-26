# Coefficient Derivatives of an MCP Penalty

[`penalty_gradient()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
returns MCP's defining derivative carried back through the map, and
[`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
its second derivative, which is \\-1/\gamma\\ inside the taper and zero
beyond it.

## Arguments

- pen:

  An
  [`McpPenalty()`](https://statmodels7.github.io/penalties7/reference/ScadPenalty.md)
  object.

- beta:

  A numeric vector of length `pen@n_coef`.

- theta:

  A list holding `lambda` and `gamma`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`penalty_gradient()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
a numeric vector of length `pen@n_coef`;
[`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
a symmetric `pen@n_coef` by `pen@n_coef` base matrix, diagonal under the
identity map and negative semidefinite always.

## Details

On each coordinate, with \\u = \lvert t\rvert\\ and \\s\\ its sign,

\$\$\frac{\partial\rho}{\partial t} = s\\(\lambda - u/\gamma)\_{+},
\qquad \frac{\partial^2\rho}{\partial t^2} = \begin{cases} -1/\gamma & u
\le \gamma\lambda \\ 0 & u \> \gamma\lambda. \end{cases}\$\$

The Hessian is **negative** wherever the penalty is doing anything,
which is the non-convexity MCP is named for, and it is constant there:
MCP is the penalty with the least curvature among those that threshold
at rate \\\lambda\\ and reach zero bias by \\\gamma\lambda\\. A
penalized objective stays convex only while the likelihood's own
curvature exceeds \\1/\gamma\\.

At \\t = 0\\ the sign is zero, so the gradient returned is `0`: one
element of the subdifferential \\\[-\lambda, \lambda\]\\. Selection
happens through
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md).

## See also

[`penalty_value.McpPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.McpPenalty.md)
for the quantity differentiated,
[`penalty_grad_theta.McpPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.McpPenalty.md)
for the hyperparameter blocks,
[`penalty_gradient.ScadPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.ScadPenalty.md)
for the sibling family.

## Examples

``` r
pen <- mcp_penalty(n_coef = 3)
th <- list(lambda = 1, gamma = 3)

# The gradient falls linearly and is exactly zero past gamma * lambda.
penalty_gradient(pen, c(0.5, 2, 5), th)
#> [1] 0.8333333 0.3333333 0.0000000

# The curvature is a constant -1/gamma inside the taper.
diag(penalty_hessian(pen, c(0.5, 2, 5), th))
#> [1] -0.3333333 -0.3333333  0.0000000
```
