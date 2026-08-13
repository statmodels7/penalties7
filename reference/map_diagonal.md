# The Diagonal of a Map That Has One

The entries of \\D\\ where the map is diagonal, and `NULL` where there
is no map or the map is not diagonal.

## Usage

``` r
map_diagonal(pen)
```

## Arguments

- pen:

  A
  [`penalty`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

## Value

A numeric vector, or `NULL`.

## Details

A diagonal map RESCALES each coordinate on its own, and a separable
penalty under one is still separable. That is what standardization comes
to: penalizing a column divided by its own spread is penalizing
\\\rho(s_j\beta_j)\\, so the scaling never has to touch the design and
the sparsity of a block survives it. A general \\D\\ mixes coordinates
and turns the problem into the generalized-lasso one, which is a
different algorithm rather than a different formula.

The map is recognized by its CLASS and not by inspecting its entries: a
Matrix diagonal object says what it is and costs \\q\\ numbers, where
testing a dense matrix for diagonality would cost \\q^2\\ and defeat the
point.
