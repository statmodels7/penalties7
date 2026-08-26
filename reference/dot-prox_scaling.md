# The Diagonal a Proximal Operator Rescales By

Returns the map's diagonal where there is one, `NULL` where there is no
map, and an error naming the penalty where the map is not diagonal. The
three answers are the three things a proximal method must do: rescale,
proceed, or decline.

## Usage

``` r
.prox_scaling(pen)
```

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

## Value

A numeric vector of diagonal entries, or `NULL` when the map is absent.

## Details

The identity behind the rescaling is

\$\$\mathrm{prox}\_{t\rho(d\\\cdot)}(v) = \mathrm{prox}\_{t d^2 \rho}(d
v)/d,\$\$

so a caller with a diagonal map applies \\d\\ to the point, squares it
into the step, calls the identity-map formula and divides back.

The error is raised here, at the call, rather than left to a formula
that would silently ignore the map.
[`has_prox()`](https://statmodels7.github.io/penalties7/reference/has_prox.md)
answers `FALSE` for the same penalty, so a caller that asks first never
reaches it.

## Errors

A message naming the penalty and saying that a non-diagonal map makes
this the generalized-lasso problem, where the operator does not split by
coordinate.

## See also

[`map_diagonal()`](https://statmodels7.github.io/penalties7/reference/map_diagonal.md),
[`.undiag()`](https://statmodels7.github.io/penalties7/reference/dot-undiag.md),
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
