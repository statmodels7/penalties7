# The Block-Diagonal Middle Matrix of a Blockwise Parent

Assembles \\\partial^2\ell/\partial b\partial b'\\ over the blocks. The
parent returns one \\p \times p\\ matrix when its response Hessian does
not depend on the observation, as the gaussian's does not, and a \\p
\times p \times n\\ array when it does.

## Usage

``` r
dp_blockdiag(pen, h, nblk)
```

## Arguments

- pen:

  A `DistribPenalty` object.

- h:

  The parent's `distrib_hess_y`.

- nblk:

  The number of blocks.

## Value

A symmetric matrix of side `nblk * pen@block`.
