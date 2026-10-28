# Apply the Linear Map and Its Transpose

`map_apply()` computes \\t = D\beta\\, carrying a coefficient vector to
the argument at which \\\rho\\ is evaluated. `map_back()` computes
\\D'g\\, carrying a gradient in \\t\\ back to a gradient in \\\beta\\. A
`NULL` map is the identity and both return their argument untouched.

## Usage

``` r
map_apply(penalty, beta)

map_back(penalty, g)
```

## Arguments

- penalty:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object, whose `map` is \\D\\ with \\m\\ rows and \\q\\ columns, or
  `NULL`.

- beta:

  A numeric vector of length \\q\\. Only `map_apply()` takes this
  argument.

- g:

  A numeric vector of length \\m\\. Only `map_back()` takes this
  argument.

## Value

`map_apply()` returns a numeric vector of length \\m\\; `map_back()` a
numeric vector of length \\q\\. Both are plain numeric vectors even when
the map is a Matrix object, rather than one-column Matrix objects.
