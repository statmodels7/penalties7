# Carry a Diagonal Middle Matrix Through the Map

\\D' \mathrm{diag}(h) D\\ for the separable Hessians, without forming
the diagonal matrix.

## Usage

``` r
map_quad(pen, h)
```

## Arguments

- pen:

  A
  [`penalty`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- h:

  A numeric vector of diagonal entries.

## Value

A `q x q` symmetric matrix.
