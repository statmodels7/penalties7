# Assemble a Piecewise Linear Table

Builds the three matrices that
[`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
returns from a function that gives one column per piece for a single
coefficient. `prox_table()` takes a function whose only argument is the
step, for a family whose pieces are the same in every coordinate;
`prox_table_j()` takes a function of the step and the coordinate index,
for a family whose rate can differ by coordinate, as the rate
\\\lambda/c_j\\ of SCAD and MCP does.

## Usage

``` r
prox_table(step, n_coef, pieces)

prox_table_j(step, n_coef, pieces)
```

## Arguments

- step:

  The step lengths, a numeric vector recycled to `n_coef`.

- n_coef:

  How many coefficients, a single whole number.

- pieces:

  A function returning a matrix of three rows, `cut`, `slope` and
  `icept`, and one column per piece. For `prox_table()` it takes one
  argument, the step; for `prox_table_j()` it takes two, the step and
  the coordinate index \\j\\.

## Value

A list of three matrices, `cut`, `slope` and `icept`, each `n_coef` by
the piece count.

## Details

The step vector is recycled to `n_coef`, and the function is called once
per coefficient with that coefficient's own step. The builders of the
branches differ only in that function, so the recycling and the
transposition are written here once. The number of pieces is read from
the first coefficient's columns, and the others must have the same
number, which holds for every shipped builder because the piece count
depends on the family and not on the step or the coordinate.

## See also

[`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md),
[`prox_apply()`](https://statmodels7.github.io/penalties7/reference/prox_apply.md)
