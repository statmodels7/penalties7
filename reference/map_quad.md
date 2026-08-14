# Carry a Middle Matrix Through the Map

\\D' \mathrm{diag}(h) D\\ for the separable Hessians, without forming
the diagonal matrix, and \\D' M D\\ for a parent read blockwise, whose
middle matrix is block diagonal rather than diagonal.

## Usage

``` r
map_quad(pen, h)

map_quad_full(pen, m)
```

## Arguments

- pen:

  A
  [`penalty`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- h:

  A numeric vector of diagonal entries.

- m:

  A symmetric matrix.

## Value

A `q x q` symmetric matrix.
