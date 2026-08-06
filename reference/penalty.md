# S7 Base Class for Penalties

The abstract parent of every penalty. A penalty is
\\\rho(D\beta;\theta)\\: a linear map \\D\\, a scalar function \\\rho\\
and hyperparameters \\\theta\\. Each hyperparameter travels with a
linkfunctions7 link, so a consumer can optimize on the unconstrained
scale.

## Usage

``` r
penalty(
  penalty_name = character(0),
  map = NULL,
  n_coef = integer(0),
  params = character(0),
  params_bounds = list(),
  link_params = list(),
  params_smooth = logical(0)
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

## Value

An object inheriting from class `penalty`.

## See also

[`quadratic_penalty`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md),
[`distrib_penalty`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md),
[`scad_penalty`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)

## Examples

``` r
S7::S7_inherits(quadratic_penalty(diag(3)), penalty)
#> [1] TRUE
```
