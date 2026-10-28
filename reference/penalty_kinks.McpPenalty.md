# Smoothness and Kind of an MCP Penalty

[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
returns the points where MCP changes branch, and
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
returns `FALSE`, since MCP is bounded and so not a density at any
constant.

## Arguments

- penalty:

  An
  [`McpPenalty()`](https://statmodels7.github.io/penalties7/reference/ScadPenalty.md)
  object.

- theta:

  A list holding `lambda` and `gamma`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
returns a numeric vector of distinct points: `0`, then
\\-\gamma\lambda/c_j\\ and \\\gamma\lambda/c_j\\ over the distinct rates
\\\lambda/c_j\\. At unit curvature these are the three points `0`,
\\-\gamma\lambda\\, \\\gamma\lambda\\, in this order.
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
returns the single logical `FALSE`.

## Details

At unit curvature the points are `0` and \\\pm\gamma\lambda\\, and they
differ in the order of the derivative that jumps. At the origin the
first derivative jumps from \\-\lambda\\ to \\\lambda\\, which is the
kink that produces exact zeros. At \\\pm\gamma\lambda\\ the first
derivative is continuous, because both branches give zero, and the
second jumps from \\-1/\gamma\\ to \\0\\ (from \\-c_j/\gamma\\ with
curvature \\c_j\\). MCP has one pair of outer points where SCAD has two,
because its taper starts at the origin instead of at \\\lambda\\.

With `curvature`, the points \\\pm\gamma\lambda/c_j\\ are returned for
each distinct value of \\\lambda/c_j\\, so the vector has more than
three points when the curvatures differ. The outer points move with the
hyperparameters, so `theta` is required.

## See also

[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
and
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
for the generics,
[`penalty_kinks.ScadPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.ScadPenalty.md)
for the sibling family.

## Examples

``` r
penalty_kinks(mcp_penalty(), list(lambda = 1, gamma = 3))
#> [1]  0 -3  3
penalty_kinks(mcp_penalty(), list(lambda = 0.5, gamma = 3))
#> [1]  0.0 -1.5  1.5

# Two curvatures give two rates and five points.
penalty_kinks(mcp_penalty(n_coef = 2, curvature = c(1, 2)),
              list(lambda = 1, gamma = 3))
#> [1]  0.0 -3.0 -1.5  3.0  1.5
is_proper(mcp_penalty())
#> [1] FALSE
```
