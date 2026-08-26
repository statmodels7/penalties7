# The Piecewise Regions of SCAD and MCP

Split the mapped coefficients into the regions each family's formulas
are written over, and return the sign and the absolute value alongside.
Shared by every method of the two branches, so a formula reads as one
[`ifelse()`](https://rdrr.io/r/base/ifelse.html) over indicators instead
of recomputing the comparisons.

## Usage

``` r
scad_parts(t, lam, a)

mcp_parts(t, lam, gam)
```

## Arguments

- t:

  The mapped coefficients \\D\beta\\, a numeric vector.

- lam:

  The hyperparameter \\\lambda\\, a single positive number.

- a:

  SCAD's shape parameter, a single number above 2. `scad_parts()` only.

- gam:

  MCP's shape parameter, a single number above 1. `mcp_parts()` only.

## Value

`scad_parts()` a list of five vectors as long as `t`: `s` the sign, `u`
the absolute value, and the logical indicators `r1` (\\u \le \lambda\\),
`r2` (the taper) and `r3` (\\u \> a\lambda\\). `mcp_parts()` a list of
three: `s`, `u`, and `r1` for \\u \le \gamma\lambda\\, the complement
needing no name.

## Details

SCAD has three regions and MCP has two, which is the whole difference
between the two helpers. Every derivative of either family is even in
\\t\\ up to the sign carried out front, so the formulas are written in
`u = abs(t)` and multiplied by `s = sign(t)` where an odd order needs
it.

At `t = 0` the sign is `0`, so a coefficient sitting exactly at the kink
gets a gradient of zero: one element of the subdifferential
\\\[-\lambda, \lambda\]\\, chosen for being the one a smooth method can
use.

## See also

[`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md),
[`penalty_value.ScadPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.ScadPenalty.md),
[`penalty_value.McpPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.McpPenalty.md)
