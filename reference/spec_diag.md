# Carry a Table Through a Diagonal Map

Builds the table of a separable penalty under a diagonal map from the
builder of its identity-map table, and returns `NULL` where the penalty
has a map that is not recognized as diagonal. With no map the builder is
called directly.

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

  A numeric vector of step lengths, of length one or `pen@n_coef`,
  recycled to `pen@n_coef`.

- build:

  A function of a penalty and a step returning the table, or `NULL`. It
  is called on the penalty with its map removed.

## Value

The list that
[`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
returns, or `NULL` when the map is not recognized as diagonal or the
builder returns `NULL`.

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

The step passed to the builder is \\t d^2\\, and the step limit of SCAD
and MCP is tested on it, so a standardized penalty takes shorter steps
and the table can be `NULL` at a step where the table without the map
exists.

The map is recognized by its class through
[`map_diagonal()`](https://statmodels7.github.io/penalties7/reference/map_diagonal.md),
which returns `NULL` for a base matrix and for a diagonal with a zero
entry. The constructors
[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md),
[`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
and
[`mcp_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
store a diagonal base matrix as a
[`Matrix::Diagonal()`](https://rdrr.io/pkg/Matrix/man/Diagonal.html)
object, so its diagonal map is recognized.

## See also

[`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md),
[`map_diagonal()`](https://statmodels7.github.io/penalties7/reference/map_diagonal.md),
[`prox_table()`](https://statmodels7.github.io/penalties7/reference/prox_table.md)
