# Coefficient Derivatives of an MCP Penalty

[`penalty_gradient()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
returns MCP's defining derivative carried back through the map, and
[`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
returns its second derivative carried back through the map. At unit
curvature the second derivative is \\-1/\gamma\\ inside the taper and
zero beyond it.

## Arguments

- penalty:

  An
  [`McpPenalty()`](https://statmodels7.github.io/penalties7/reference/ScadPenalty.md)
  object.

- beta:

  A numeric vector of length `penalty@n_coef`.

- theta:

  A list holding `lambda` and `gamma`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`penalty_gradient()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
returns a numeric vector of length `penalty@n_coef`;
[`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
a symmetric `penalty@n_coef` by `penalty@n_coef` base matrix, diagonal
under the identity map and negative semidefinite always.

## Details

With unit curvature, on each coordinate, with \\u = \lvert t\rvert\\ and
\\s\\ its sign,

\$\$\frac{\partial\rho}{\partial t} = s\\(\lambda - u/\gamma)\_{+},
\qquad \frac{\partial^2\rho}{\partial t^2} = \begin{cases} -1/\gamma & u
\le \gamma\lambda \\ 0 & u \> \gamma\lambda. \end{cases}\$\$

With curvature \\c_j\\ the taper ends at \\\gamma\lambda/c_j\\, the
first derivative is \\s\\(\lambda - c_j u/\gamma)\_{+}\\, and the second
derivative inside the taper is \\-c_j/\gamma\\.

The Hessian is negative wherever the penalty is not constant, which is
the non-convexity in the name of MCP, and it is constant there. Among
the penalties whose derivative is \\\lambda\\ at the origin and zero
beyond \\\gamma\lambda\\, MCP has the smallest maximum concavity, the
supremum of \\-\rho''\\ (Zhang 2010). A penalized objective stays convex
only while the curvature of the likelihood in each coordinate exceeds
\\1/\gamma\\ (\\c_j/\gamma\\ with `curvature`).

At \\t = 0\\ the sign is zero, so the gradient returned is `0`, one
element of the subdifferential \\\[-\lambda, \lambda\]\\. Exact zeros
are produced by
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
