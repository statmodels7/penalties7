# The Tolerance of the Rank Identity

Returns `1e-3`, the relative departure from
\\\sum_k\lambda_k\operatorname{tr}(S^{+}P_k) = r\\ beyond which
[`additive_sum()`](https://statmodels7.github.io/penalties7/reference/additive_sum.md)
reports `NaN` instead of a number. The departure is taken relative to
the rank \\r\\.

## Usage

``` r
additive_check_tol()
```

## Value

A single number.

## Details

The identity is exact, so a departure from it is accumulated rounding in
the similarity transformation of
[`additive_stable()`](https://statmodels7.github.io/penalties7/reference/additive_stable.md).
The threshold is chosen above the departures that rounding produces at
points that the transformation resolves, so that a departure beyond it
signals that the smaller eigenvalues have been lost. The threshold is
not a bound on the error of the returned quantities: at a point that
passes the check, the error of the log pseudo-determinant can exceed the
departure from the identity.

## See also

[`additive_sum()`](https://statmodels7.github.io/penalties7/reference/additive_sum.md),
which reads it.
