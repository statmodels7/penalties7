# S7 Classes for the Derivative-Defined Penalties

The classes
[`scad_penalty`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
and
[`mcp_penalty`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
instantiate: families the literature defines by \\\rho'\\ and whose
\\\rho\\ follows by the closed piecewise antiderivative anchored at
\\\rho(0) = 0\\. Both are improper: \\\rho\\ is bounded, so
\\\exp(-\rho)\\ does not integrate.

## Usage

``` r
ScadPenalty(
  penalty_name = character(0),
  map = NULL,
  n_coef = integer(0),
  params = character(0),
  params_bounds = list(),
  link_params = list(),
  params_smooth = logical(0)
)

McpPenalty(
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

An object of class `ScadPenalty` or `McpPenalty`.

## See also

[`scad_penalty`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md),
[`mcp_penalty`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)

## Examples

``` r
S7::S7_inherits(scad_penalty(), ScadPenalty)
#> [1] TRUE
```
