# Align and Validate the Hyperparameters

Puts a hyperparameter argument into the single shape that every branch
reads: a list reordered to `pen@params`, with the names removed from the
values, and checked against `pen@params_bounds` as open intervals.
Returns the aligned list. For a penalty with no hyperparameters it
returns an empty list without examining `theta`.

## Usage

``` r
align_ptheta(pen, theta)
```

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- theta:

  A named list of hyperparameter values, or a named numeric vector
  carrying the same. Extra entries are dropped; a missing one is an
  error naming which. Each value must be a single number.

## Value

A list of the same length and order as `pen@params`, each element an
unnamed single number. [`list()`](https://rdrr.io/r/base/list.html) when
the penalty has no hyperparameters.

## Details

A named numeric vector is accepted in place of a list and converted with
[`as.list()`](https://rdrr.io/r/base/list.html), because the SCAD and
MCP methods read their hyperparameters with `$`, which does not apply to
an atomic vector. The conversion is made here, in the function that
every generic calls, so every branch receives a list.

The bounds are **open**, so a hyperparameter at an endpoint is rejected
rather than clamped: `alpha = 1` on an elastic net whose bound is \\(0,
1)\\ signals an error. Parameters in distributions7 are validated in the
same way, and the open bounds keep the value of every link finite.

## Errors

`Missing parameter(s) in 'theta': ... Expected: ...` when a name is
absent or `theta` is unnamed, `Parameter 'p' must be a single number.`
when a value is `NULL`, empty or longer than one, and
`Parameter 'p' must lie in the open interval (a, b).` when a value is
non-finite or outside its bounds.
