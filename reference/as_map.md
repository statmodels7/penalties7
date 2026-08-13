# A Map, in Whatever Form It Keeps

The map as the caller gave it, densified only where it is not already a
matrix of some kind.

## Usage

``` r
as_map(map)
```

## Arguments

- map:

  A matrix, a Matrix, or anything coercible to one.

## Value

The map.

## Details

A Matrix object is kept as it is. Densifying a diagonal map would cost
\\q^2\\ numbers where it holds \\q\\, and a diagonal map is exactly what
standardization is: a rescaling of each coordinate, under which a
separable penalty stays separable and its proximal operator stays
closed. Every arithmetic the map takes part in – the product, the
crossproduct – works for both kinds.
