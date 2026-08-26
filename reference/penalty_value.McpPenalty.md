# Value of an MCP Penalty

Returns the closed antiderivative of MCP's defining \\\rho'\\, summed
over the coordinates of \\D\beta\\ and anchored at \\\rho(0) = 0\\. No
constant is added, MCP being improper.

## Arguments

- pen:

  An
  [`McpPenalty()`](https://statmodels7.github.io/penalties7/reference/ScadPenalty.md)
  object.

- beta:

  A numeric vector of length `pen@n_coef`, already coerced by the
  generic.

- theta:

  A list holding `lambda` and `gamma`, already aligned and bound-checked
  by the generic.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A single number, the sum over coordinates. Bounded above by
`pen@n_coef * gamma * lambda^2 / 2`.

## Details

With \\u = \lvert t\rvert\\ on each coordinate,

\$\$\rho(u) = \lambda u - \frac{u^2}{2\gamma} \\\\ (u \le
\gamma\lambda), \qquad \rho(u) = \frac{\gamma\lambda^2}{2} \\\\ (u \>
\gamma\lambda),\$\$

the two pieces meeting continuously at \\\gamma\lambda\\. Where SCAD
holds the lasso's slope over an interval before tapering, MCP starts
tapering at once, so its value is a downward parabola from the origin.

## See also

[`mcp_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
for the construction,
[`penalty_gradient.McpPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.McpPenalty.md)
for the derivative it integrates,
[`penalty_value.ScadPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.ScadPenalty.md)
for the sibling family.

## Examples

``` r
pen <- mcp_penalty(n_coef = 1)
th <- list(lambda = 1, gamma = 3)

# Inside the taper and beyond it.
sapply(c(0.5, 2, 100), function(t) penalty_value(pen, t, th))
#> [1] 0.4583333 1.3333333 1.5000000

# The flat value is gamma lambda^2 / 2.
3 / 2
#> [1] 1.5
```
