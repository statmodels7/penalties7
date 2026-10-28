# A Penalty With Its Map Removed

Returns the same penalty with `map` set to `NULL`, which is the object
that the diagonal-map transport passes to the formulas written for the
identity map.

## Usage

``` r
.undiag(penalty)
```

## Arguments

- penalty:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object whose map is diagonal.

## Value

The same object with `map` set to `NULL`.

## Details

The transport applies the diagonal to the point and to the step before
the call and divides the result by the diagonal afterwards, so the
formula needs the penalty without its map. The other properties are
unchanged; in particular `n_coef` is unaffected, because a diagonal map
is square.

## See also

[`.prox_scaling()`](https://statmodels7.github.io/penalties7/reference/dot-prox_scaling.md),
[`spec_diag()`](https://statmodels7.github.io/penalties7/reference/spec_diag.md),
[`map_diagonal()`](https://statmodels7.github.io/penalties7/reference/map_diagonal.md)
