# Value of an MCP Penalty

Returns the closed antiderivative of MCP's defining \\\rho'\\, summed
over the coordinates of \\D\beta\\ and anchored at \\\rho(0) = 0\\. It
adds no constant, since MCP is improper.

## Arguments

- penalty:

  An
  [`McpPenalty()`](https://statmodels7.github.io/penalties7/reference/ScadPenalty.md)
  object.

- beta:

  A numeric vector of length `penalty@n_coef`, already coerced by the
  generic.

- theta:

  A list holding `lambda` and `gamma`, already aligned and bound-checked
  by the generic.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A single number, the sum over coordinates. It is bounded above by
\\\sum_j \gamma\lambda^2/(2c_j)\\, where the sum runs over the
coordinates of \\D\beta\\; under the identity map at unit curvature the
bound is `penalty@n_coef * gamma * lambda^2 / 2`.

## Details

With unit curvature and \\u = \lvert t\rvert\\ on each coordinate,

\$\$\rho(u) = \lambda u - \frac{u^2}{2\gamma} \\\\ (u \le
\gamma\lambda), \qquad \rho(u) = \frac{\gamma\lambda^2}{2} \\\\ (u \>
\gamma\lambda),\$\$

where the two pieces meet continuously at \\\gamma\lambda\\. SCAD keeps
the slope of the lasso over an interval before it tapers, while MCP
tapers from the origin, so its value is a downward parabola up to
\\\gamma\lambda\\. With curvature \\c_j\\ (the `curvature` property),
coordinate \\j\\ contributes \\c_j\\\rho(u_j;\\ \lambda/c_j)\\, the same
formula with \\\lambda/c_j\\ in place of \\\lambda\\ and multiplied by
\\c_j\\.

## See also

[`mcp_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
for the construction,
[`penalty_gradient.McpPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.McpPenalty.md)
for the derivative that it integrates,
[`penalty_value.ScadPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.ScadPenalty.md)
for the sibling family.

## Examples

``` r
pen <- mcp_penalty(n_coef = 1)
th <- list(lambda = 1, gamma = 3)

# The coordinate is evaluated twice inside the taper and once beyond it.
sapply(c(0.5, 2, 100), function(t) penalty_value(pen, t, th))
#> [1] 0.4583333 1.3333333 1.5000000

# The flat value is gamma lambda^2 / 2.
3 / 2
#> [1] 1.5
```
