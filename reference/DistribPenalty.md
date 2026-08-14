# S7 Class for the Separable Penalty

The class
[`distrib_penalty`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md)
instantiates: a penalty built by applying a univariate distributions7
log-density coordinatewise to \\D\beta\\.

## Usage

``` r
DistribPenalty(
  penalty_name = character(0),
  map = NULL,
  n_coef = integer(0),
  params = character(0),
  params_bounds = list(),
  link_params = list(),
  params_smooth = logical(0),
  parent = NULL,
  kinks = integer(0)
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

- parent:

  The distributions7 object.

- kinks:

  The declared non-differentiable points of the parent's log-density in
  its argument.

## Value

An object of class `DistribPenalty`.

## See also

[`distrib_penalty`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md)

## Examples

``` r
S7::S7_inherits(ridge_penalty(n_coef = 2), DistribPenalty)
#> [1] FALSE
```
