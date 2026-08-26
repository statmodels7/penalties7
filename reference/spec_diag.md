# Carry a Table Through a Diagonal Map

Builds the table of a separable penalty under a diagonal map from the
builder of its identity-map table, and returns `NULL` where there is a
map and it is not diagonal. With no map at all the builder is called
directly.

## Usage

``` r
spec_diag(pen, step, build)
```

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object whose map is `NULL` or diagonal.

- step:

  A numeric vector of step lengths, recycled to `pen@n_coef`.

- build:

  A function of a penalty and a step returning the table, or `NULL`.
  Called on the penalty with its map removed.

## Value

The list
[`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
returns, or `NULL` when the map is not diagonal or the builder itself
declines.

## Details

A diagonal map only rescales each coordinate, and the identity

\$\$\mathrm{prox}\_{t\rho(d\\\cdot)}(v) = \mathrm{prox}\_{t d^2 \rho}(d
v)/d\$\$

carries the table across. The identity-map table reads \\\lvert w\rvert
\le \mathrm{cut}\\ to \\\mathrm{sign}(w)(\mathrm{slope}\\\lvert
w\rvert + \mathrm{icept})\\ at \\w = dv\\, so a cut on \\\lvert
w\rvert\\ is a cut on \\\lvert v\rvert\\ divided by \\\lvert d\rvert\\,
the intercept divides by the same, and the slope does not move,
multiplying a point that was scaled and then divided back. The operator
is odd, so only the magnitude of \\d\\ enters.

The step handed to the builder is \\t d^2\\, which is where the
convexity condition of SCAD and MCP tightens: a standardized penalty
takes shorter steps and a table may come back `NULL` where the unmapped
one exists.

The map is recognized by its class through
[`map_diagonal()`](https://statmodels7.github.io/penalties7/reference/map_diagonal.md),
so a base matrix that happens to be diagonal is not treated as one.

## See also

[`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md),
[`map_diagonal()`](https://statmodels7.github.io/penalties7/reference/map_diagonal.md),
[`prox_table()`](https://statmodels7.github.io/penalties7/reference/prox_table.md)
