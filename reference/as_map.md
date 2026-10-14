# A Map, in Whatever Storage It Arrived In

Returns the map unchanged when it is already a Matrix object, and
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) of it otherwise.
Called once, by each constructor that takes a map, so that a map given
as a data frame or a vector becomes a matrix while a sparse or diagonal
one keeps its own storage.

## Usage

``` r
as_map(map)
```

## Arguments

- map:

  A matrix, a Matrix, or anything
  [`as.matrix()`](https://rdrr.io/r/base/matrix.html) accepts.

## Value

The same object when it is a Matrix, and a base matrix otherwise.

## Details

A dense copy of a diagonal map would hold \\q^2\\ numbers instead of
\\q\\. A diagonal map is what standardization produces, a rescaling of
each coordinate under which a separable penalty stays separable and its
proximal operator stays in closed form. The product and the
crossproduct, which are the only operations applied to the map, are
defined for both kinds of storage.
