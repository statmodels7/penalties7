# The Diagonal of a Map That Has One

Returns the diagonal entries of \\D\\ where the map is a Matrix diagonal
object, and `NULL` where there is no map, where the map is of another
class, or where a diagonal entry is zero or missing.

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

A diagonal map rescales each coordinate on its own, and a separable
penalty under one is still separable. That is what standardization comes
to: penalizing a column divided by its own spread is penalizing
\\\rho(s_j\beta_j)\\, so the scaling never has to touch the design and
the sparsity of a block survives it. A general \\D\\ mixes coordinates
and turns the problem into the generalized-lasso one, which is a
different algorithm and not a different formula.

The map is recognized by its **class** and not by inspecting its
entries: a Matrix diagonal object says what it is and costs \\q\\
numbers, where testing a dense matrix for diagonality would cost \\q^2\\
and defeat the point. A base matrix that happens to be diagonal
therefore answers `NULL`, and a penalty carrying one has no proximal
operator.

A zero on the diagonal is rejected as well, the transport dividing by
it.

## See also

[`has_prox()`](https://statmodels7.github.io/penalties7/reference/has_prox.md),
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md),
[`spec_diag()`](https://statmodels7.github.io/penalties7/reference/spec_diag.md)
