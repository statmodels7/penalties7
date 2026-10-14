# S7 Class for the Structured Quadratic Penalty

The class of the objects that
[`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md)
builds: the Gaussian prior whose precision is a parameters7 matrix
parameter. It adds one property to
[`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md),
the structure itself, and every quantity of the branch is computed from
the methods of the structure.

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
  [`print()`](https://rdrr.io/r/base/print.html) and by any code that
  reports which penalty a block carries.

- map:

  The matrix \\D\\, of \\m\\ rows and \\q\\ columns, or `NULL` for the
  identity. A Matrix object is kept in its own storage, so a diagonal
  map, which is what standardization produces, is recognized as diagonal
  by its class.

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

  A parameters7 `matrix_parameter`, read as the prior's precision. It
  supplies the dimension, the rank, the null basis, the free names, the
  matrix and its derivative arrays.

## Value

An S7 object of class `StructuredPenalty`, inheriting from
[`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md),
with the seven inherited properties and `structure`.

## Details

The hyperparameters are the structure's free vector, so `params` is its
`free_names`, every entry of `params_bounds` is \\(-\infty, \infty)\\
and every link is the identity. The constraint that keeps the matrix
positive definite is part of the structure, where a scalar link cannot
express it. `map` is always `NULL`.

## Methods

The methods registered on this class, by the page that documents each.

- The value and its coefficient derivatives:
  [penalty_value()](https://statmodels7.github.io/penalties7/reference/penalty_value.StructuredPenalty.md),
  [penalty_gradient()](https://statmodels7.github.io/penalties7/reference/penalty_gradient.StructuredPenalty.md)
  and
  [penalty_hessian()](https://statmodels7.github.io/penalties7/reference/penalty_gradient.StructuredPenalty.md).

- The hyperparameter derivatives:
  [penalty_grad_theta()](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.StructuredPenalty.md),
  [penalty_hess_theta()](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.StructuredPenalty.md)
  and
  [penalty_cross()](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.StructuredPenalty.md).

- The predicates:
  [penalty_kinks()](https://statmodels7.github.io/penalties7/reference/penalty_kinks.StructuredPenalty.md),
  [is_proper()](https://statmodels7.github.io/penalties7/reference/penalty_kinks.StructuredPenalty.md)
  and
  [is_quadratic()](https://statmodels7.github.io/penalties7/reference/penalty_kinks.StructuredPenalty.md).

- The quantities of a marginal criterion:
  [penalty_matrix()](https://statmodels7.github.io/penalties7/reference/penalty_matrix.StructuredPenalty.md),
  [penalty_rank()](https://statmodels7.github.io/penalties7/reference/penalty_matrix.StructuredPenalty.md),
  [penalty_null_basis()](https://statmodels7.github.io/penalties7/reference/penalty_matrix.StructuredPenalty.md)
  and
  [penalty_logpdet()](https://statmodels7.github.io/penalties7/reference/penalty_matrix.StructuredPenalty.md).

- The higher derivatives of a marginal criterion:
  [penalty_dhessian()](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.StructuredPenalty.md),
  [penalty_d2hessian()](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.StructuredPenalty.md),
  [penalty_dcross()](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.StructuredPenalty.md)
  and
  [beta_quadratic()](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.StructuredPenalty.md).

- The movement of the Hessian with the coefficients:
  [penalty_dhessian_beta()](https://statmodels7.github.io/penalties7/reference/penalty_dhessian_beta.penalty.md),
  [penalty_d2hessian_beta()](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.penalty.md)
  and
  [penalty_dhessian_beta_theta()](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.penalty.md).
  The last two return zero matrices, the Hessian being free of the
  coefficients, and are documented on the page of the base-class
  methods.

- The proximal operator:
  [penalty_prox()](https://statmodels7.github.io/penalties7/reference/penalty_prox.QuadraticPenalty.md),
  documented on the page of the quadratic branch.

- The interpretable quantities:
  [penalty_readable()](https://statmodels7.github.io/penalties7/reference/penalty_readable.StructuredPenalty.md),
  which returns those that the structure declares through
  [`parameters7::param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.html).

The class registers no
[`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
or
[`penalty_draw()`](https://statmodels7.github.io/penalties7/reference/penalty_draw.md)
method, so the methods of the base class
[`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
answer those generics.

## See also

[`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md)
for the constructor to use,
[`penalty_value.StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.StructuredPenalty.md)
for what the branch computes,
[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
for the branch whose matrix is fixed.

## Examples

``` r
pen <- structured_penalty(parameters7::ar1(4))
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
