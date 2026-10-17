# The Quadratic Form of a Quadratic Penalty

Computes \\(D\beta)'P(D\beta)\\, the part of the value that depends on
the coefficients. It is shared by
[`penalty_value.QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.QuadraticPenalty.md)
and the hyperparameter derivatives, which use \\\beta\\ only through
this quantity.

## Usage

``` r
quad_form(pen, beta)
```

## Arguments

- pen:

  A
  [`QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/QuadraticPenalty.md)
  object.

- beta:

  A numeric vector of length `pen@n_coef`.

## Value

A single number, non-negative for a positive semidefinite \\P\\.

## Details

Uses \\P\\ and the map, not the cached \\D'PD\\, so the intermediate is
the \\m\\-vector \\D\beta\\ and the function does not form a \\q \times
q\\ product.

## See also

[`penalty_value.QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.QuadraticPenalty.md)
