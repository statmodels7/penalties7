# The Diagonal of a Map That Has One

Returns the diagonal entries of \\D\\ when the map is a Matrix diagonal
object, and `NULL` when there is no map, when the map has another class,
or when a diagonal entry is zero or missing.

## Usage

``` r
map_diagonal(pen)
```

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

## Value

A numeric vector of length `nrow(pen@map)`, or `NULL`. A unit diagonal
object (`diag = "U"`) gives a vector of ones.

## Details

A diagonal map rescales each coordinate on its own, so a separable
penalty under a diagonal map stays separable. Standardization is such a
map: penalizing a column divided by its own spread \\s_j\\ is penalizing
\\\rho(s_j\beta_j)\\, so the scaling leaves the design unchanged and the
sparsity of a block is kept. A general \\D\\ mixes the coordinates and
turns the problem into the generalized-lasso problem, which requires a
different algorithm.

The map is recognized by its class, without inspecting its entries: a
Matrix diagonal object stores \\q\\ numbers, while testing a dense
matrix for diagonality would cost \\q^2\\ operations at every call. A
base matrix that is diagonal in value therefore gives `NULL`.
[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md),
[`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
and
[`mcp_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
convert a square map whose off-diagonal entries are all zero to a
[`Matrix::Diagonal()`](https://rdrr.io/pkg/Matrix/man/Diagonal.html)
when the penalty is built, so a separable, SCAD or MCP penalty built
from such a base matrix has a proximal operator.

A zero on the diagonal also gives `NULL`, because the transport divides
by the diagonal entries.

## See also

[`has_prox()`](https://statmodels7.github.io/penalties7/reference/has_prox.md),
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md),
[`spec_diag()`](https://statmodels7.github.io/penalties7/reference/spec_diag.md)
