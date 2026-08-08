# S7 Class for a Sum of Quadratic Penalties

The class
[`additive_penalty`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
instantiates: several quadratic penalties added together, each with a
smoothing parameter of its own.

## Usage

``` r
AdditivePenalty(
  penalty_name = character(0),
  map = NULL,
  n_coef = integer(0),
  params = character(0),
  params_bounds = list(),
  link_params = list(),
  params_smooth = logical(0),
  mats = list(),
  p_rank = integer(0)
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

- mats:

  The list of penalty matrices, already carried through the map.

- p_rank:

  The rank of the sum, the same at every positive parameter value.

## Value

An object of class `AdditivePenalty`.

## See also

[`additive_penalty`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)

## Examples

``` r
S7::S7_inherits(additive_penalty(list(diag(3), diag(c(1, 0, 0)))),
                AdditivePenalty)
#> [1] TRUE
```
