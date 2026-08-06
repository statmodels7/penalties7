# S7 Class for the Quadratic Penalty

The class
[`quadratic_penalty`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
instantiates. Beyond the base properties it stores the matrix \\P\\, its
rank, its null basis and its log pseudo-determinant, all fixed by one
eigendecomposition at construction.

## Usage

``` r
QuadraticPenalty(
  penalty_name = character(0),
  map = NULL,
  n_coef = integer(0),
  params = character(0),
  params_bounds = list(),
  link_params = list(),
  params_smooth = logical(0),
  P = NULL,
  p_rank = integer(0),
  null_basis = NULL,
  logpdet_P = integer(0),
  DPD = NULL
)
```

## Arguments

- penalty_name:

  A string naming the penalty.

- map:

  The matrix \\D\\, or `NULL` for the identity.

- n_coef:

  The number of coefficients \\q\\.

- params:

  Hyperparameter names, in order.

- params_bounds:

  A named list of open intervals.

- link_params:

  A named list of linkfunctions7 links.

- params_smooth:

  Logical vector; which hyperparameters are differentiable.

- P:

  The symmetric positive semidefinite matrix.

- p_rank:

  The rank of \\P\\.

- null_basis:

  An orthonormal basis of the null space of \\D'PD\\.

- logpdet_P:

  The log pseudo-determinant of \\P\\.

- DPD:

  The assembled \\D'PD\\, cached.

## Value

An object of class `QuadraticPenalty`.

## See also

[`quadratic_penalty`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)

## Examples

``` r
S7::S7_inherits(quadratic_penalty(diag(2)), QuadraticPenalty)
#> [1] TRUE
```
