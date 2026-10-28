# The Zero Matrix of an Additive Penalty's Side

Returns the zero matrix of side `penalty@n_coef`: sparse when the
penalty is repeated blockwise, so that a third-order block added to it
by a caller is not allocated dense, and a base matrix otherwise.

## Usage

``` r
additive_zero(penalty)
```

## Arguments

- penalty:

  An
  [`AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/AdditivePenalty.md)
  object.

## Value

A square zero matrix.
