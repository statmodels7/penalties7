# The Diagonal by Which a Proximal Operator Rescales

Returns the diagonal of the map when the map is a Matrix diagonal object
with nonzero entries, `NULL` when there is no map, and signals an error
naming the penalty for any other map. The three outcomes correspond to
the three cases of a proximal method: the point is rescaled, the
operator is applied as it is, or the call stops.

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

so a caller with a diagonal map multiplies the point by \\d\\ and the
step by \\d^2\\, calls the formula written for the identity map, and
divides the result by \\d\\.

The error is signaled here because the formulas written for the identity
map would otherwise ignore the map.
[`has_prox()`](https://statmodels7.github.io/penalties7/reference/has_prox.md)
returns `FALSE` for the same penalty, so a caller that tests it first
does not reach the error.

## Errors

The message names the penalty and states that its map is not a Matrix
diagonal object with nonzero entries, and that the operator is available
only under the identity map and such a diagonal map, because under a
general map it does not split by coordinate.

## See also

[`map_diagonal()`](https://statmodels7.github.io/penalties7/reference/map_diagonal.md),
[`.undiag()`](https://statmodels7.github.io/penalties7/reference/dot-undiag.md),
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
