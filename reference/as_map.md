# A Map, in Whatever Storage It Arrived In

Returns the map unchanged when it is already a Matrix object, and
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) of it otherwise.
Called once, by each branch's constructor, so that a map given as a data
frame or a vector becomes a matrix while a sparse or diagonal one keeps
its own storage.

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

Densifying a diagonal map would cost \\q^2\\ numbers where it holds
\\q\\, and a diagonal map is exactly what standardization is: a
rescaling of each coordinate, under which a separable penalty stays
separable and its proximal operator stays closed. Every arithmetic the
map takes part in, the product and the crossproduct, is defined for both
kinds.
