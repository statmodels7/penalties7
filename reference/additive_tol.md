# The Split Threshold of the Components

Returns \\\epsilon^{0.4}\\ (about `5.5e-07`), with \\\epsilon\\ the
machine precision. In
[`additive_stable()`](https://statmodels7.github.io/penalties7/reference/additive_stable.md),
a component whose contribution
\\\lambda_k\max\_{ij}\lvert(P_k)\_{ij}\rvert\\ is at least this fraction
of the largest contribution joins the dominant group, and the other
components are carried to the next level.

## Usage

``` r
additive_tol()
```

## Value

A single number.

## Details

The threshold balances two sources of rounding. A smaller value widens
the dominant group, so one decomposition must resolve a spread of up to
the reciprocal of the threshold. A larger value narrows the group, so
the recursion runs deeper and every level adds the rounding of one more
reduction. The exponent 0.4 was chosen on the rank identity of
[`additive_check_tol()`](https://statmodels7.github.io/penalties7/reference/additive_check_tol.md);
mgcv's `gam.reparam` uses \\\epsilon^{0.3}\\ in the same role.

## See also

[`additive_stable()`](https://statmodels7.github.io/penalties7/reference/additive_stable.md),
which reads it.
