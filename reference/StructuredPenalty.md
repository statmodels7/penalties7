# S7 Class for the Structured Quadratic Penalty

The class
[`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md)
builds: the Gaussian prior whose covariance or precision is a
parameters7 matrix parameter. It adds one property to
[`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md),
the structure itself, and every quantity the branch supplies is read off
that structure's own contract.

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

  A single string naming the penalty, used by
  [`print()`](https://rdrr.io/r/base/print.html) and by consumers that
  report which penalty a block carries.

- map:

  The matrix \\D\\, of \\m\\ rows and \\q\\ columns, or `NULL` for the
  identity. A Matrix object is kept in its own storage; a diagonal map
  is what standardization comes to and is recognized by its class.

- n_coef:

  The number of coefficients \\q\\. A single whole number.

- params:

  The hyperparameter names, in the order every derivative list is keyed
  by. `character(0)` for a penalty with none.

- params_bounds:

  A named list, one entry per hyperparameter, each a numeric pair giving
  an **open** interval. A value at either endpoint is rejected, so
  `(0, Inf)` excludes zero.

- link_params:

  A named list, one linkfunctions7 link per hyperparameter, carrying
  that hyperparameter's own interval onto the whole real line.

- params_smooth:

  A logical vector, one entry per hyperparameter, `TRUE` where the value
  is differentiable in it.

- structure:

  A parameters7 `matrix_parameter` whose `role` is `"covariance"` or
  `"precision"`. It supplies the dimension, the rank, the null basis,
  the free names, the matrix and its derivative arrays.

## Value

An S7 object of class `StructuredPenalty`, inheriting from
[`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md),
with the seven inherited properties and `structure`.

## Details

The hyperparameters are the structure's free vector, so `params` is its
`free_names`, every entry of `params_bounds` is \\(-\infty, \infty)\\
and every link is the identity. The constraint that keeps the matrix
positive definite lives inside the structure, where a scalar link cannot
express it. `map` is always `NULL`.

## See also

[`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md)
for the constructor to use,
[`penalty_value.StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.StructuredPenalty.md)
for what the branch computes,
[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
for the branch whose matrix is fixed.

## Examples

``` r
pen <- structured_penalty(parameters7::ar1(4, role = "precision"))
S7::S7_inherits(pen, StructuredPenalty)
#> [1] TRUE

# The hyperparameters are the structure's free names, unconstrained and
# on identity links.
pen@params
#> [1] "log_scale" "z_rho"    
pen@params_bounds
#> $log_scale
#> [1] -Inf  Inf
#> 
#> $z_rho
#> [1] -Inf  Inf
#> 
```
