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

  A numeric vector of length \\m\\, the diagonal entries. `map_quad()`
  only.

- m:

  A symmetric \\m \times m\\ matrix. `map_quad_full()` only.

## Value

A \\q \times q\\ symmetric base matrix, from both.

## Details

Both return a base matrix even when the map is a Matrix. A Matrix map
carries its class through the crossproduct, and the result would then be
the one thing in the contract that is not a base matrix: the
identity-map branch is already dense at any width,
[`map_back()`](https://statmodels7.github.io/penalties7/reference/map_apply.md)
coerces its vector for the same reason, and a consumer writing this into
a block of its own information fails on the class before it fails on the
arithmetic. The coercion is done here, where the contract is stated.
