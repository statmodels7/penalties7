# The Diagonal and the Log-Jacobian of a Map

`dp_diag_entries()` returns the entries of a Matrix diagonal map,
including the zeros, and `NULL` for any other map. `map_log_jacobian()`
returns \\\log\lvert\det D\rvert\\: `0` for a `NULL` map,
\\\sum_j\log\lvert d_j\rvert\\ for a diagonal one, the log-modulus of
the determinant for a square nonsingular one, and `NA` for a map that is
not square or is singular.
[`dp_log_jacobian()`](https://statmodels7.github.io/penalties7/reference/dp_log_jacobian.md)
returns the same quantity for a built penalty.

## Usage

``` r
dp_diag_entries(map)

map_log_jacobian(map)
```

## Arguments

- map:

  The map \\D\\: `NULL`, a base matrix or a Matrix.

## Value

`dp_diag_entries()` returns a numeric vector of length `nrow(map)`, a
unit diagonal (`diag = "U"`) giving ones, or `NULL`.
`map_log_jacobian()` returns a single number, possibly `NA`.

## Details

A diagonal map is treated as singular when an entry is zero or missing.
Any other square map is converted once to a dense matrix; it is treated
as singular when an entry is not finite or when the rank that
[`qr()`](https://rdrr.io/r/base/qr.html) reports at its default
tolerance is below the number of columns, and otherwise the log-modulus
comes from [`determinant()`](https://rdrr.io/r/base/det.html).

## See also

[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md),
[`penalty_value.DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.DistribPenalty.md)
