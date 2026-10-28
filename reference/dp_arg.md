# The Parent's Argument, Shaped and Unshaped

`dp_arg()` reshapes \\D\beta\\ into the argument of the parent: the
vector itself for a univariate parent, and a matrix of one row per block
for a multivariate one. `dp_flat()` undoes the reshaping, carrying the
result of the parent back to a vector in coefficient order.

## Usage

``` r
dp_arg(penalty, t)

dp_flat(penalty, g)
```

## Arguments

- penalty:

  A
  [`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
  object.

- t:

  The mapped coefficient vector \\D\beta\\, of length a multiple of
  `penalty@block`. This argument belongs to `dp_arg()` only.

- g:

  The result of the parent: a vector for a univariate parent, a matrix
  of one row per block for a multivariate one. This argument belongs to
  `dp_flat()` only.

## Value

`dp_arg()` returns the vector unchanged when `penalty@block` is `1L`,
and otherwise a matrix with `penalty@block` columns filled by row.
`dp_flat()` returns a numeric vector, in coefficient order.

## Details

The blocks are the successive stretches of \\D\beta\\, so the reshaping
fills by row: block \\i\\ occupies positions \\(i-1)p+1, \dots, ip\\.
This is the order in which a grouped design assembles its coefficients
and the order that \\I_m \otimes \Sigma\\ assumes, so the two agree.

Both functions return their argument unchanged when the block width is
one, which is the case of every univariate parent.

## See also

[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md),
[`dp_blockdiag()`](https://statmodels7.github.io/penalties7/reference/dp_blockdiag.md)
for the matrix-valued counterpart.
