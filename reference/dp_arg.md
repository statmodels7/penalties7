# The Parent's Argument, Shaped and Unshaped

`dp_arg` reshapes \\D\beta\\ into the argument the parent reads – the
vector itself for a univariate parent, one row per block for a
multivariate one – and `dp_flat` undoes it.

## Usage

``` r
dp_arg(pen, t)

dp_flat(pen, g)
```

## Arguments

- pen:

  A `DistribPenalty` object.

- t:

  The mapped coefficient vector.

- g:

  The parent's answer, a vector or a matrix of one row per block.

## Value

A vector or a matrix.

## Details

The blocks are the SUCCESSIVE stretches of \\D\beta\\, so the reshaping
fills by row: block \\i\\ occupies positions \\(i-1)p+1, \dots, ip\\.
That is the order a grouped design assembles its coefficients in, and
the order \\I_m \otimes \Sigma\\ assumes.
