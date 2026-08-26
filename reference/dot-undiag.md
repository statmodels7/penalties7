# A Penalty With Its Map Removed

Returns the same penalty with `map` set to `NULL`, which is the object
the diagonal-map transport hands to the identity-map formulas.

## Usage

``` r
.undiag(pen)
```

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object whose map is diagonal.

## Value

The same object with `map` set to `NULL`.

## Details

The transport applies the diagonal to the point and to the step before
calling, and divides the answer back afterwards, so what the formula
needs is the penalty without its map. Nothing else on the object
changes, and in particular `n_coef` is unaffected, a diagonal map being
square.

## See also

[`.prox_scaling()`](https://statmodels7.github.io/penalties7/reference/dot-prox_scaling.md),
[`spec_diag()`](https://statmodels7.github.io/penalties7/reference/spec_diag.md),
[`map_diagonal()`](https://statmodels7.github.io/penalties7/reference/map_diagonal.md)
