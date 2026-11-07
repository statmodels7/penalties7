# The Log-Jacobian of a Separable Penalty's Map

Returns \\\log\lvert\det D\rvert\\ for the map \\D\\ of a separable
penalty: `0` when the penalty has no map, the sum of the logarithms of
the absolute diagonal entries for a diagonal map, and the stored
`log_jacobian` property for any other map.

## Usage

``` r
dp_log_jacobian(penalty)
```

## Arguments

- penalty:

  A
  [`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
  object.

## Value

A single number, possibly `NA`.

## Details

For a diagonal map the quantity is recomputed from the entries, at the
cost of \\q\\ logarithms, so that the result stays correct when a caller
replaces the map with a subset of its diagonal through
[`S7::set_props()`](https://rconsortium.github.io/S7/reference/props.html)
and the stored property is out of date. A diagonal map with a zero or
missing entry gives `NA`, as does a stored property that is `NA` because
the map is not square or is singular.

## See also

[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md)

## Examples

``` r
d <- distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0)
pen <- distrib_penalty(d, n_coef = 3)
dp_log_jacobian(pen)
#> [1] 0
pen2 <- distrib_penalty(d, map = Matrix::Diagonal(x = c(1, 2, 4)))
all.equal(dp_log_jacobian(pen2), log(8))
#> [1] TRUE
```
