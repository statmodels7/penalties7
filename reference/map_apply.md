# Apply the Linear Map and Its Transpose

`map_apply()` computes \\t = D\beta\\, carrying a coefficient vector to
the argument \\\rho\\ is evaluated at. `map_back()` computes \\D'g\\,
carrying a gradient in \\t\\ back to a gradient in \\\beta\\. A `NULL`
map is the identity and both return their argument untouched.

## Usage

``` r
map_apply(pen, beta)

map_back(pen, g)
```

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object, whose `map` is \\D\\ with \\m\\ rows and \\q\\ columns, or
  `NULL`.

- beta:

  A numeric vector of length \\q\\. `map_apply()` only.

- g:

  A numeric vector of length \\m\\. `map_back()` only.

## Value

`map_apply()` a numeric vector of length \\m\\; `map_back()` a numeric
vector of length \\q\\. Both are plain numeric even when the map is a
Matrix, so a consumer never meets a one-column `Matrix` where it
expected a vector.
