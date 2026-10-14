# S7 Class for a Sum of Quadratic Penalties

The class of the objects that
[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
builds: several quadratic penalties added together, each with a
smoothing parameter of its own. Beyond the seven properties that every
[`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
carries, an object of this class stores the list of component matrices,
already carried through the map, the rank and null basis of the sum, and
the eigenvalue tolerance at which ranks are counted.

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
  p_rank = integer(0),
  null_basis = NULL,
  tol = 1e-10
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

- mats:

  A list of symmetric positive semidefinite matrices of the same side,
  already carried through the map, so each is `n_coef` by `n_coef`. One
  smoothing parameter goes with each, named `lambda1`, `lambda2` and so
  on in this order.

- p_rank:

  The rank of the sum, a single whole number, the same at every positive
  setting of the parameters.

- null_basis:

  An orthonormal basis of the null space shared by the components,
  `n_coef` by `n_coef - p_rank`, read from the same eigendecomposition
  as `p_rank` and fixed in the same way.

- tol:

  The relative eigenvalue tolerance, a single positive number, `1e-10`
  by default.
  [`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
  counts `p_rank` at this tolerance, and
  [`additive_sum()`](https://statmodels7.github.io/penalties7/reference/additive_sum.md)
  passes it to
  [`additive_stable()`](https://statmodels7.github.io/penalties7/reference/additive_stable.md)
  as `r_tol`, so the rank of each group of components is counted at the
  same tolerance.

## Value

An S7 object of class `AdditivePenalty`, inheriting from
[`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md),
with the seven inherited properties and the four above. Its `map` is
always `NULL`, a map given to the constructor having been absorbed into
`mats`.

## Details

The rank and the null basis are stored because both are properties of
the components alone. The null space of a sum of positive semidefinite
matrices with positive weights is the intersection of the null spaces of
the components, so neither the rank nor the null basis changes with the
smoothing parameters. A rank counted on the assembled \\S(\lambda)\\
would fall as the parameters spread apart, because the smaller
contributions drop below any fixed relative eigenvalue tolerance; the
last example of
[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
shows this.

## Methods

The methods registered on this class, by the page that documents each.

- The value and its coefficient derivatives:
  [penalty_value()](https://statmodels7.github.io/penalties7/reference/penalty_value.AdditivePenalty.md),
  [penalty_gradient()](https://statmodels7.github.io/penalties7/reference/penalty_gradient.AdditivePenalty.md)
  and
  [penalty_hessian()](https://statmodels7.github.io/penalties7/reference/penalty_gradient.AdditivePenalty.md).

- The hyperparameter derivatives:
  [penalty_grad_theta()](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.AdditivePenalty.md),
  [penalty_hess_theta()](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.AdditivePenalty.md)
  and
  [penalty_cross()](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.AdditivePenalty.md).

- The predicates:
  [penalty_kinks()](https://statmodels7.github.io/penalties7/reference/penalty_kinks.AdditivePenalty.md),
  [is_proper()](https://statmodels7.github.io/penalties7/reference/penalty_kinks.AdditivePenalty.md)
  and
  [is_quadratic()](https://statmodels7.github.io/penalties7/reference/is_quadratic.AdditivePenalty.md).

- The quantities of a marginal criterion:
  [penalty_matrix()](https://statmodels7.github.io/penalties7/reference/penalty_matrix.AdditivePenalty.md),
  [penalty_rank()](https://statmodels7.github.io/penalties7/reference/penalty_matrix.AdditivePenalty.md),
  [penalty_null_basis()](https://statmodels7.github.io/penalties7/reference/penalty_matrix.AdditivePenalty.md)
  and
  [penalty_logpdet()](https://statmodels7.github.io/penalties7/reference/penalty_matrix.AdditivePenalty.md).

- The higher derivatives of a marginal criterion:
  [penalty_dhessian()](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.AdditivePenalty.md),
  [penalty_d2hessian()](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.AdditivePenalty.md),
  [penalty_dcross()](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.AdditivePenalty.md)
  and
  [beta_quadratic()](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.AdditivePenalty.md).

- The movement of the Hessian with the coefficients:
  [penalty_dhessian_beta()](https://statmodels7.github.io/penalties7/reference/penalty_dhessian_beta.penalty.md),
  [penalty_d2hessian_beta()](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.penalty.md)
  and
  [penalty_dhessian_beta_theta()](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.penalty.md).
  The last two return zero matrices, the Hessian being free of the
  coefficients, and are documented on the page of the base-class
  methods.

The class registers no
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md),
[`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md),
[`penalty_draw()`](https://statmodels7.github.io/penalties7/reference/penalty_draw.md)
or
[`penalty_readable()`](https://statmodels7.github.io/penalties7/reference/penalty_readable.md)
method, so calls to those generics dispatch to the methods of the base
class
[`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md).

## See also

[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
for the constructor to use,
[`penalty_value.AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.AdditivePenalty.md)
for the value of the penalty,
[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
for one matrix with one scale.

## Examples

``` r
pen <- additive_penalty(list(diag(3), diag(c(1, 0, 0))))
S7::S7_inherits(pen, AdditivePenalty)
#> [1] TRUE

# One hyperparameter per component, and a rank fixed at construction.
pen@params
#> [1] "lambda1" "lambda2"
c(rank = pen@p_rank, null_width = ncol(pen@null_basis))
#>       rank null_width 
#>          3          0 
length(pen@mats)
#> [1] 2
```
