# S7 Class for the Structured Quadratic Penalty

The class
[`structured_penalty`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md)
instantiates: the Gaussian prior whose covariance or precision is a
parameters7 matrix parameter, so that the hyperparameters enter the
matrix itself.

## Usage

``` r
StructuredPenalty(
  penalty_name = character(0),
  map = NULL,
  n_coef = integer(0),
  params = character(0),
  params_bounds = list(),
  link_params = list(),
  params_smooth = logical(0),
  structure = NULL
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

- structure:

  The parameters7 `matrix_parameter`.

## Value

An object of class `StructuredPenalty`.

## See also

[`structured_penalty`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md)

## Examples

``` r
S7::S7_inherits(
  structured_penalty(parameters7::log_cholesky(2, role = "precision")),
  StructuredPenalty)
#> [1] TRUE
```
