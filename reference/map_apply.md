# Apply the Linear Map and Its Transpose

`map_apply` computes \\t = D\beta\\ and `map_back` computes \\D'g\\; a
`NULL` map is the identity and pays nothing.

## Usage

``` r
map_apply(pen, beta)

map_back(pen, g)
```

## Arguments

- pen:

  A
  [`penalty`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- beta:

  A numeric vector of coefficients.

- g:

  A numeric vector of length `nrow(D)`.

## Value

A numeric vector.
