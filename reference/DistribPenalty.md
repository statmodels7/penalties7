# S7 Class for the Separable Penalty

The class of the penalties that
[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md)
builds: a penalty formed by applying a distributions7 log-density to the
successive blocks of \\D\beta\\. Beyond the seven properties that every
penalty carries, the class stores the parent distribution, the points
where the parent's log-density has a kink, the block width (one for a
univariate parent and the dimension for a multivariate one) and the
log-Jacobian of the map.

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
  kinks = integer(0),
  block = integer(0),
  log_jacobian = NA_real_
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

- parent:

  The distributions7 object that supplies the log-density, usually a
  [`distributions7::fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.html)
  wrapper holding the location at zero.

- kinks:

  A numeric vector of the points where the parent's log-density is not
  differentiable in its argument, possibly empty.
  [`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md)
  derives it with
  [`distrib_kinks()`](https://statmodels7.github.io/penalties7/reference/distrib_kinks.md)
  unless its argument `kinks` is supplied.

- block:

  The block width: `1L` for a univariate parent and the parent's `n_dim`
  for a multivariate one.

- log_jacobian:

  \\\log\lvert\det D\rvert\\, computed by
  [map_log_jacobian()](https://statmodels7.github.io/penalties7/reference/dp_diag_entries.md)
  when the penalty is built: `0` without a map, a finite number for a
  square nonsingular map, and `NA` for a map that is not square or is
  singular. The default is `NA`.

## Value

An S7 object of class `DistribPenalty`, inheriting from
[`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md),
with the seven inherited properties and the four above.

## Details

The hyperparameters are the parent's own free parameters. The properties
`params`, `params_bounds`, `link_params` and `params_smooth` are copied
from the distribution object, so the hyperparameters of the penalty,
their bounds, their links and their smoothness are those of the parent.

## Methods

The methods registered on this class, by the page that documents each.

- The value and its coefficient derivatives:
  [penalty_value()](https://statmodels7.github.io/penalties7/reference/penalty_value.DistribPenalty.md),
  [penalty_gradient()](https://statmodels7.github.io/penalties7/reference/penalty_gradient.DistribPenalty.md)
  and
  [penalty_hessian()](https://statmodels7.github.io/penalties7/reference/penalty_gradient.DistribPenalty.md).

- The hyperparameter derivatives:
  [penalty_grad_theta()](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.DistribPenalty.md),
  [penalty_hess_theta()](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.DistribPenalty.md)
  and
  [penalty_cross()](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.DistribPenalty.md).

- The predicates and the structure of the prior:
  [penalty_kinks()](https://statmodels7.github.io/penalties7/reference/penalty_kinks.DistribPenalty.md),
  [is_proper()](https://statmodels7.github.io/penalties7/reference/penalty_kinks.DistribPenalty.md),
  [penalty_null_basis()](https://statmodels7.github.io/penalties7/reference/penalty_null_basis.DistribPenalty.md)
  and
  [penalty_readable()](https://statmodels7.github.io/penalties7/reference/penalty_readable.DistribPenalty.md).

- The higher derivatives of a marginal criterion:
  [penalty_dhessian()](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.DistribPenalty.md),
  [penalty_d2hessian()](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.DistribPenalty.md),
  [penalty_dcross()](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.DistribPenalty.md)
  and
  [beta_quadratic()](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.DistribPenalty.md).

- The movement of the Hessian with the coefficients:
  [penalty_dhessian_beta()](https://statmodels7.github.io/penalties7/reference/penalty_dhessian_beta.penalty.md),
  [penalty_d2hessian_beta()](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.penalty.md)
  and
  [penalty_dhessian_beta_theta()](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.penalty.md).

- The proximal operator and the draws:
  [penalty_prox()](https://statmodels7.github.io/penalties7/reference/penalty_prox.DistribPenalty.md),
  [penalty_prox_spec()](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
  and
  [penalty_draw()](https://statmodels7.github.io/penalties7/reference/penalty_draw.DistribPenalty.md).

The class registers no
[is_quadratic()](https://statmodels7.github.io/penalties7/reference/is_quadratic.md),
[penalty_matrix()](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md),
[penalty_rank()](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
or
[penalty_logpdet()](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
method, so the methods of the base class
[`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
are used for those generics.

## See also

[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md)
for the constructor,
[`ridge_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md)
and its siblings for the named instances,
[`penalty_value.DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.DistribPenalty.md)
for the value of the penalty.

## Examples

``` r
pen <- lasso_penalty(n_coef = 2)
S7::S7_inherits(pen, DistribPenalty)
#> [1] TRUE

# The parent is a fixed Laplace, the kink is at zero, the blocks are single
# coordinates, and the hyperparameter is the rate of the parent.
pen@parent
#> Distribution: Fixed Laplace2 [mu=0]
#> Type:         Continuous
#> Dimensions:   univariate
#> 
#> Parameters:
#>   lambda (rate)               | Link: log        | Domain: (0, Inf)
#> 
#> Fixed:
#>   mu = 0
pen@kinks
#> [1] 0
pen@block
#> [1] 1
pen@params
#> [1] "lambda"

# Without a map the log-Jacobian is zero.
pen@log_jacobian
#> [1] 0
```
