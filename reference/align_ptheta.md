# Align and Validate the Hyperparameters

Puts a hyperparameter argument into the single shape that every branch
reads: a list reordered to `penalty@params`, whose values are unnamed,
and checked against `penalty@params_bounds` as open intervals. Returns
the aligned list. For a penalty with no hyperparameters it returns an
empty list without examining `theta`.

## Usage

``` r
align_ptheta(penalty, theta)
```

## Arguments

- penalty:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- theta:

  A named list of hyperparameter values, or a named numeric vector
  carrying the same. Extra entries are dropped; a missing entry is an
  error that names the entry. Each value must be a single number.

## Value

A list of the same length and order as `penalty@params`, in which each
element is an unnamed single number. The value is
[`list()`](https://rdrr.io/r/base/list.html) when the penalty has no
hyperparameters.

## Details

A named numeric vector is accepted in place of a list and converted with
[`as.list()`](https://rdrr.io/r/base/list.html), because the SCAD and
MCP methods read their hyperparameters with `$`, which does not apply to
an atomic vector. This function, which every generic calls, makes the
conversion, so every branch receives a list.

The bounds are open, so a hyperparameter at an endpoint is rejected
rather than clamped: `alpha = 1` on an elastic net whose bound is \\(0,
1)\\ signals an error. Parameters in distributions7 are validated in the
same way, and the open bounds keep the value of every link finite.

## Errors

The error message is
`'theta' lacks the parameter(s) ...; the expected parameters are ....`
when a name is absent or `theta` is unnamed,
`Parameter 'p' must be a single number.` when a value is `NULL`, empty
or longer than one, and
`Parameter 'p' must lie in the open interval (a, b).` when a value is
non-finite or outside its bounds.
