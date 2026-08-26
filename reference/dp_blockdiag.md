# The Block-Diagonal Middle Matrix of a Blockwise Parent

Assembles the parent's \\\partial^2\ell/\partial b\partial b'\\ into the
block-diagonal matrix a multivariate separable penalty's Hessian is
sandwiched with: one \\p \times p\\ block per block of coefficients, and
zeros between them.

## Usage

``` r
dp_blockdiag(pen, h, nblk)
```

## Arguments

- pen:

  A
  [`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
  object with `pen@block` above one.

- h:

  The parent's
  [`distributions7::distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.html):
  a \\p \times p\\ matrix or a \\p \times p \times n\\ array.

- nblk:

  The number of blocks, a single whole number.

## Value

A symmetric base matrix of side `nblk * pen@block`, block diagonal with
`nblk` blocks.

## Details

The zeros between the blocks are the separability. Coordinates in
different blocks are independent under the prior, so the Hessian has no
entry linking them, however dependent the coordinates within a block
are.

Two input shapes are accepted because distributions7 returns two. A
family whose response Hessian does not depend on the observation, as the
multivariate gaussian's does not, returns one \\p \times p\\ matrix and
it is placed in every block; a family whose does returns a \\p \times p
\times n\\ array and each slice goes to its own block. The two are told
apart with [`is.matrix()`](https://rdrr.io/r/base/matrix.html).

## See also

[`dp_arg()`](https://statmodels7.github.io/penalties7/reference/dp_arg.md),
and
[`map_quad_full()`](https://statmodels7.github.io/penalties7/reference/map_quad.md),
which carries this through the map.
