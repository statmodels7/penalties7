# The Parent's Argument, Shaped and Unshaped

`dp_arg()` reshapes \\D\beta\\ into the argument the parent reads: the
vector itself for a univariate parent, and a matrix of one row per block
for a multivariate one. `dp_flat()` undoes the reshaping, carrying the
parent's answer back to a vector in coefficient order.

## Usage

``` r
dp_arg(pen, t)

dp_flat(pen, g)
```

## Arguments

- pen:

  A
  [`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
  object.

- t:

  The mapped coefficient vector \\D\beta\\, of length a multiple of
  `pen@block`. `dp_arg()` only.

- g:

  The parent's answer: a vector for a univariate parent, a matrix of one
  row per block for a multivariate one. `dp_flat()` only.

## Value

`dp_arg()` the vector unchanged when `pen@block` is `1L`, and otherwise
a matrix with `pen@block` columns filled by row. `dp_flat()` a numeric
vector, in coefficient order.

## Details

The blocks are the **successive** stretches of \\D\beta\\, so the
reshaping fills by row: block \\i\\ occupies positions \\(i-1)p+1,
\dots, ip\\. That is the order a grouped design assembles its
coefficients in, and the order \\I_m \otimes \Sigma\\ assumes, so the
two agree without either being told about the other.

Both are the identity when the block width is one, which is every
univariate parent, and the pair costs nothing there.

## See also

[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md),
[`dp_blockdiag()`](https://statmodels7.github.io/penalties7/reference/dp_blockdiag.md)
for the matrix-valued counterpart.
