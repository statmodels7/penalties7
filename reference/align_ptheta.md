# Align and Validate the Hyperparameters

Puts a hyperparameter argument into the one shape every branch reads:
reordered to `pen@params`, with stray names stripped off the values, and
checked against `pen@params_bounds` treated as open intervals. Returns
the aligned list. A penalty with no hyperparameters returns an empty
list without looking at `theta`.

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
  error naming which. Each value may be a vector, in which case every
  element is bound-checked.

## Value

A list of the same length and order as `pen@params`, each element
unnamed. [`list()`](https://rdrr.io/r/base/list.html) when the penalty
has no hyperparameters.

## Details

A named numeric vector carries what the list carries, and the branches
split on how they read it: `[[` accepts both, `$` accepts only the list.
A caller passing a vector therefore reached the quadratic and separable
branches and failed inside SCAD and MCP, three frames down and naming
neither the argument nor the penalty. Converting here, at the one point
every generic passes through, settles the shape for all of them.

The bounds are **open**, so a hyperparameter at an endpoint is rejected
rather than clamped: `alpha = 1` on an elastic net whose bound is \\(0,
1)\\ throws. That matches distributions7, whose parameters are validated
the same way, and it is what keeps a link's inverse finite.

## Errors

`Missing parameter(s) in 'theta': ... Expected: ...` when a name is
absent or `theta` is unnamed, and
`Parameter 'p' must lie in the open interval (a, b).` when a value is
non-finite or outside its bounds.
