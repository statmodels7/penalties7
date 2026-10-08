# The Curvature of Each Coordinate of a SCAD or MCP Penalty

Validates the curvature a SCAD or MCP penalty is scaled by, and returns
it one value per coordinate of \\D\beta\\.

## Usage

``` r
check_curv(curv, n)

curv_of(pen, n)
```

## Arguments

- curv:

  `NULL`, or positive finite numbers, one per coordinate of \\D\beta\\
  or a single value repeated.

- n:

  The number of coordinates of \\D\beta\\.

- pen:

  A
  [`ScadPenalty()`](https://statmodels7.github.io/penalties7/reference/ScadPenalty.md)
  or
  [`McpPenalty()`](https://statmodels7.github.io/penalties7/reference/ScadPenalty.md).

## Value

`check_curv()` a numeric vector of length `n`, or `numeric(0)` for
`NULL`; `curv_of()` a numeric vector of length `n`.

## Details

SCAD and MCP are defined on the canonical problem \\\tfrac12(z -
\theta)^2 + p\_\lambda(\lvert\theta\rvert)\\, whose loss has unit
curvature. Where the loss has curvature \\c_j\\ in coordinate \\j\\, the
penalty is written \\c_j\\p(t_j;\\ \lambda/c_j,\\ a)\\: the slope at the
origin stays \\\lambda\\, and the knee moves to \\a\lambda/c_j\\, so
that the ratio of the knee to the soft threshold is \\a\\ in every
coordinate. `curv_of()` reads the stored curvature, and gives ones where
none is stored, which is the penalty as defined without this scaling.
