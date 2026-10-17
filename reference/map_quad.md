# Carry a Middle Matrix Through the Map

`map_quad()` computes \\D' \mathrm{diag}(h) D\\ without forming the
diagonal matrix, which is the Hessian of a separable penalty carried
back to the coefficients. `map_quad_full()` computes \\D'MD\\ for a
parent read blockwise, whose middle matrix is block diagonal. With a
`NULL` map the first returns `diag(h)` and the second returns \\M\\.

## Usage

``` r
map_quad(pen, h)

map_quad_full(pen, m)
```

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object, whose `map` is \\D\\ with \\m\\ rows and \\q\\ columns, or
  `NULL`.

- h:

  A numeric vector of length \\m\\, the diagonal entries. Only
  `map_quad()` takes this argument.

- m:

  A symmetric \\m \times m\\ matrix. Only `map_quad_full()` takes this
  argument.

## Value

A \\q \times q\\ symmetric base matrix, from both.

## Details

Both return a base matrix even when the map is a Matrix object. The
crossproduct of a Matrix map would otherwise be a Matrix object, whereas
the result for a `NULL` map is a base matrix and
[`map_back()`](https://statmodels7.github.io/penalties7/reference/map_apply.md)
returns a plain vector. `map_quad()` and `map_quad_full()` convert the
result, so that the caller receives a base matrix from every branch and
can write it into a block of a base information matrix.
