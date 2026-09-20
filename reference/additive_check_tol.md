# How Far the Rank Identity May Be Violated

The relative departure from \\\sum_k\lambda_k\operatorname{tr}(S^{+}P_k)
= r\\ past which
[`additive_sum()`](https://statmodels7.github.io/penalties7/reference/additive_sum.md)
reports `NaN` rather than a number.

## Usage

``` r
additive_check_tol()
```

## Value

A single number.

## Details

It is a threshold on an exact identity and not an accuracy claim.
Measured over 200 random parameter vectors per shape on nine shapes –
anisotropic `te()` at two and three margins and with margins of
different dimension, and `basis7::adaptive_smooth()` at three, five,
eight and twelve components – with each parameter drawn log-uniformly
over the range a fit visits, the worst violation at
[`additive_tol()`](https://statmodels7.github.io/penalties7/reference/additive_tol.md)
is `5.3e-05` and the median is at machine precision. The threshold is
two orders above that worst case, so it fires where the transformation
has genuinely run out of precision and not where it has merely lost its
last digits.

## See also

[`additive_sum()`](https://statmodels7.github.io/penalties7/reference/additive_sum.md),
which reads it.
