# The Curvature of Each Coordinate of a SCAD or MCP Penalty

`check_curv()` validates the curvature by which a SCAD or MCP penalty is
scaled and returns it with one value per coordinate of \\D\beta\\.
`curv_of()` returns the curvature stored in a penalty.

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

`check_curv()` returns a numeric vector of length `n`, or `numeric(0)`
for `NULL`; `curv_of()` the stored numeric vector, or a vector of `n`
ones.

## Details

SCAD and MCP are defined on the canonical problem \\\tfrac12(z -
\theta)^2 + \rho(\theta)\\, whose loss has unit curvature. Where the
loss has curvature \\c_j\\ in coordinate \\j\\, the penalty of that
coordinate is \\c_j\\\rho(t_j;\\ \lambda/c_j)\\: the slope at the origin
stays \\\lambda\\, and the knee moves to \\a\lambda/c_j\\
(\\\gamma\lambda/c_j\\ for MCP), so that the ratio of the knee to
\\\lambda/c_j\\ is the shape parameter in every coordinate.

`check_curv()` signals an error for an empty vector, a missing,
infinite, zero or negative value, and a length other than one or `n`; a
single value is repeated `n` times. `curv_of()` returns the stored
curvature, or `n` ones where the curvature is not stored, which is the
penalty as defined without this scaling. A stored curvature is returned
whatever `n` is.
