# S7 Classes for the Derivative-Defined Penalties

`ScadPenalty` is the class
[`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
builds and `McpPenalty` the one
[`mcp_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
builds. Both are families the literature defines by \\\rho'\\ rather
than by \\\rho\\, and whose value follows from the closed piecewise
antiderivative anchored at \\\rho(0) = 0\\. Neither adds a property to
[`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md):
everything the branch needs is the two hyperparameters and the map, and
every quantity is computed from them.

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

## Value

An S7 object of class `ScadPenalty` or `McpPenalty`, inheriting from
[`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
and carrying its seven properties and no others.

## Details

Both are improper. \\\rho'\\ is exactly zero beyond \\a\lambda\\ for
SCAD and beyond \\\gamma\lambda\\ for MCP, so \\\rho\\ is bounded and
flat far from the origin and \\\exp(-\rho)\\ does not integrate over
\\\mathbb{R}^q\\ at any constant.
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
returns `FALSE` for both, and the four quantities a marginal criterion
reads are unavailable: their hyperparameters have to be chosen along a
path or by cross-validation.

The two differ in which shape parameter they carry and where it is
bounded: SCAD's is `a`, on \\(2, \infty)\\, and MCP's is `gamma`, on
\\(1, \infty)\\.

## See also

[`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
and
[`mcp_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
for the constructors,
[`penalty_value.ScadPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.ScadPenalty.md)
and
[`penalty_value.McpPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.McpPenalty.md)
for what each computes,
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
for the operator a fit steps with.

## Examples

``` r
S7::S7_inherits(scad_penalty(), ScadPenalty)
#> [1] TRUE
S7::S7_inherits(mcp_penalty(), McpPenalty)
#> [1] TRUE

# Neither adds a property to the base class.
setdiff(names(S7::props(scad_penalty())),
        names(S7::props(quadratic_penalty(diag(1)))))
#> character(0)
```
