# The Piecewise Regions of SCAD and MCP

Split the mapped coefficients into the regions over which the formulas
of each family are written, and return the sign and the absolute value
with them. The value, the coefficient derivatives and the hyperparameter
derivatives of the two classes use these helpers, so each formula is one
[`ifelse()`](https://rdrr.io/r/base/ifelse.html) over indicators and the
comparisons are made once.

## Usage

``` r
scad_parts(t, lam, a)

mcp_parts(t, lam, gam)
```

## Arguments

- t:

  The mapped coefficients \\D\beta\\, a numeric vector.

- lam:

  The rate of each coordinate, \\\lambda/c_j\\ with \\c_j\\ the
  curvature of the coordinate: a positive numeric vector as long as `t`,
  or a single positive number when every curvature is one.

- a:

  SCAD's shape parameter, a single number above 2. `scad_parts()` only.

- gam:

  MCP's shape parameter, a single number above 1. `mcp_parts()` only.

## Value

`scad_parts()` a list of five vectors as long as `t`: `s` the sign, `u`
the absolute value, and the logical indicators `r1` (\\u \le\\ `lam`),
`r2` (the taper, `lam` \\\< u \le a\\ `lam`) and `r3` (\\u \> a\\
`lam`). `mcp_parts()` a list of three vectors: `s`, `u`, and `r1` for
\\u \le \gamma\\ `lam`, whose complement is the flat region.

## Details

SCAD has three regions and MCP has two, and this is the only difference
between the two helpers. Every derivative of either family is even in
\\t\\ up to the sign in front, so the formulas are written in
`u = abs(t)` and multiplied by `s = sign(t)` where an odd order needs
it.

At `t = 0` the sign is `0`, so a coefficient exactly at the kink gets a
gradient of zero. This is one element of the subdifferential
\\\[-\lambda, \lambda\]\\, and the one that a smooth method can use.

## See also

[`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md),
[`penalty_value.ScadPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.ScadPenalty.md),
[`penalty_value.McpPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.McpPenalty.md)
