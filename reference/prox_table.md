# Assemble a Piecewise Linear Table

Builds the three matrices
[`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
returns from a function that gives one column per piece at a single
step, recycling a step that is the same for every coefficient.

## Usage

``` r
prox_table(step, n_coef, pieces)
```

## Arguments

- step:

  The step lengths, a numeric vector recycled to `n_coef`.

- n_coef:

  How many coefficients, a single whole number.

- pieces:

  A function of one step returning a matrix of three rows, `cut`,
  `slope` and `icept`, and one column per piece.

## Value

A list of three matrices, `cut`, `slope` and `icept`, each `n_coef` by
the piece count.

## Details

Every branch's builder differs only in that per-step function, so this
holds the recycling and the transposition once. The number of pieces is
read from the first coefficient's answer and every other is required to
match, which is automatic: the piece count depends on the family and on
the hyperparameters, never on the step.

## See also

[`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md),
[`prox_apply()`](https://statmodels7.github.io/penalties7/reference/prox_apply.md)
