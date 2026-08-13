# Carry a Table Through a Diagonal Map

Builds the table of a separable penalty under a diagonal map from the
builder of its identity-map table, and returns `NULL` where the map is
not diagonal.

## Usage

``` r
spec_diag(pen, step, build)
```

## Arguments

- pen:

  A
  [`penalty`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- step:

  A numeric vector of step lengths.

- build:

  A function of a penalty and a step returning the table, or `NULL`.

## Value

The list
[`penalty_prox_spec`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
returns, or `NULL`.

## Details

A diagonal map only rescales each coordinate, and the identity
\$\$\mathrm{prox}\_{t\rho(d\\\cdot)}(v) = \mathrm{prox}\_{t d^2 \rho}(d
v)/d\$\$ carries the table across. The identity-map table reads \\\|w\|
\le \mathrm{cut}\\ to \\\mathrm{sign}(w)(\mathrm{slope}\\\|w\| +
\mathrm{icept})\\ at \\w = dv\\, so a cut on \\\|w\|\\ is a cut on
\\\|v\|\\ divided by \\\|d\|\\ and the intercept divides by the same,
while the slope, which multiplies a point that was scaled and then
divided back, does not move. The operator is odd, so only the magnitude
of \\d\\ enters.

The step handed to the builder is \\t d^2\\, which is where the
convexity condition of SCAD and MCP tightens: a standardized penalty
takes shorter steps.
