# Smoothness and Kind of an MCP Penalty

[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
returns the three points where MCP changes branch, and
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
returns `FALSE`, MCP being bounded and so not a density at any constant.

## Arguments

- pen:

  An
  [`McpPenalty()`](https://statmodels7.github.io/penalties7/reference/ScadPenalty.md)
  object.

- theta:

  A list holding `lambda` and `gamma`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
a numeric vector of three points, in the order `0`, \\-\gamma\lambda\\,
\\\gamma\lambda\\.
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
the single logical `FALSE`.

## Details

The three points are `0` and \\\pm\gamma\lambda\\, and they differ in
order. At the origin the first derivative jumps from \\-\lambda\\ to
\\\lambda\\, the kink that produces exact zeros. At \\\pm\gamma\lambda\\
the first derivative is continuous, both branches giving zero, and the
second jumps from \\-1/\gamma\\ to \\0\\. Both outer points move with
the hyperparameters, so `theta` is not optional.

MCP has one pair of outer points where SCAD has two, its taper starting
at the origin instead of at \\\lambda\\.

## See also

[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
and
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
for the generics,
[`penalty_kinks.ScadPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.ScadPenalty.md)
for the sibling family's five points.

## Examples

``` r
penalty_kinks(mcp_penalty(), list(lambda = 1, gamma = 3))
#> [1]  0 -3  3
penalty_kinks(mcp_penalty(), list(lambda = 0.5, gamma = 3))
#> [1]  0.0 -1.5  1.5
is_proper(mcp_penalty())
#> [1] FALSE
```
