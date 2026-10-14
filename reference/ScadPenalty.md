# S7 Classes for the SCAD and MCP Penalties

`ScadPenalty` is the class that
[`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
builds, and `McpPenalty` is the class that
[`mcp_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
builds. The literature defines both families by the derivative
\\\rho'\\, and their value is the closed piecewise antiderivative
anchored at \\\rho(0) = 0\\. Each class adds one property to
[`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md),
`curv`, the curvature of the loss in each coordinate of \\D\beta\\, and
every quantity is computed from the two hyperparameters, the map and
`curv`.

## Usage

``` r
ScadPenalty(
  penalty_name = character(0),
  map = NULL,
  n_coef = integer(0),
  params = character(0),
  params_bounds = list(),
  link_params = list(),
  params_smooth = logical(0),
  curv = integer(0)
)

McpPenalty(
  penalty_name = character(0),
  map = NULL,
  n_coef = integer(0),
  params = character(0),
  params_bounds = list(),
  link_params = list(),
  params_smooth = logical(0),
  curv = integer(0)
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

- curv:

  The curvature of the loss in each coordinate of \\D\beta\\, a numeric
  vector of positive finite values with one entry per coordinate, or
  `numeric(0)`, which stores no curvature and gives every coordinate
  unit curvature. The constructors
  [`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
  and
  [`mcp_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
  build it with
  [`check_curv()`](https://statmodels7.github.io/penalties7/reference/check_curv.md).

## Value

An S7 object of class `ScadPenalty` or `McpPenalty`, inheriting from
[`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
and carrying its seven properties and `curv`.

## Details

Both penalties are improper. The derivative \\\rho'\\ is exactly zero
beyond \\a\lambda\\ for SCAD and beyond \\\gamma\lambda\\ for MCP, so
\\\rho\\ is bounded and flat far from the origin and \\\exp(-\rho)\\
does not integrate over \\\mathbb{R}^q\\ at any constant.
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
returns `FALSE` for both, and
[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
and the three quantities beside it signal an error, so the
hyperparameters are chosen along a path or by cross-validation.

The two classes differ in the shape parameter and its bound: SCAD's is
`a`, on \\(2, \infty)\\, and MCP's is `gamma`, on \\(1, \infty)\\.

## Methods

The two classes register methods for the same generics. For
`ScadPenalty`:

- value and coefficient derivatives:
  [penalty_value()](https://statmodels7.github.io/penalties7/reference/penalty_value.ScadPenalty.md),
  [penalty_gradient()](https://statmodels7.github.io/penalties7/reference/penalty_gradient.ScadPenalty.md),
  [penalty_hessian()](https://statmodels7.github.io/penalties7/reference/penalty_gradient.ScadPenalty.md);

- hyperparameter derivatives:
  [penalty_grad_theta()](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.ScadPenalty.md),
  [penalty_hess_theta()](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.ScadPenalty.md),
  [penalty_cross()](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.ScadPenalty.md);

- points of non-smoothness and kind:
  [penalty_kinks()](https://statmodels7.github.io/penalties7/reference/penalty_kinks.ScadPenalty.md),
  [is_proper()](https://statmodels7.github.io/penalties7/reference/penalty_kinks.ScadPenalty.md);

- proximal operator:
  [penalty_prox()](https://statmodels7.github.io/penalties7/reference/penalty_prox.ScadPenalty.md),
  [penalty_prox_spec()](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md).

For `McpPenalty`:

- value and coefficient derivatives:
  [penalty_value()](https://statmodels7.github.io/penalties7/reference/penalty_value.McpPenalty.md),
  [penalty_gradient()](https://statmodels7.github.io/penalties7/reference/penalty_gradient.McpPenalty.md),
  [penalty_hessian()](https://statmodels7.github.io/penalties7/reference/penalty_gradient.McpPenalty.md);

- hyperparameter derivatives:
  [penalty_grad_theta()](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.McpPenalty.md),
  [penalty_hess_theta()](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.McpPenalty.md),
  [penalty_cross()](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.McpPenalty.md);

- points of non-smoothness and kind:
  [penalty_kinks()](https://statmodels7.github.io/penalties7/reference/penalty_kinks.McpPenalty.md),
  [is_proper()](https://statmodels7.github.io/penalties7/reference/penalty_kinks.McpPenalty.md);

- proximal operator:
  [penalty_prox()](https://statmodels7.github.io/penalties7/reference/penalty_prox.McpPenalty.md),
  [penalty_prox_spec()](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md).

Every other generic dispatches to the method of the base class
[`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md):
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
returns `FALSE`, and
[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
and the three quantities beside it signal an error.

## See also

[`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
and
[`mcp_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
for the constructors,
[`penalty_value.ScadPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.ScadPenalty.md)
and
[`penalty_value.McpPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.McpPenalty.md)
for the values,
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
for the proximal operator used in a fit.

## Examples

``` r
S7::S7_inherits(scad_penalty(), ScadPenalty)
#> [1] TRUE
S7::S7_inherits(mcp_penalty(), McpPenalty)
#> [1] TRUE

# Both add one property to the base class, the curvature per coordinate.
setdiff(names(S7::props(scad_penalty())), names(S7::props(lasso_penalty())))
#> [1] "curv"
```
