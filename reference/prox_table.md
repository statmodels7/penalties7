# Assemble a Piecewise Linear Table

Builds the three matrices from one row per piece, recycling a step that
is the same for every coefficient.

## Usage

``` r
prox_table(step, n_coef, pieces)
```

## Arguments

- step:

  The step lengths.

- n_coef:

  How many coefficients.

- pieces:

  A function of one step returning a matrix whose rows are `cut`,
  `slope`, `icept`.

## Value

The list
[`penalty_prox_spec`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
returns.
