# The Gap at Which the Components Are Split

A component whose contribution \\\lambda_k\lVert P_k\rVert\\ is at least
this fraction of the largest joins the dominant group; the rest are
carried to the next level.

## Usage

``` r
additive_tol()
```

## Value

A single number.

## Details

The constant is measured here and not taken from elsewhere, and the
sweep that chose it has a clear interior optimum. A smaller value widens
the dominant group, so one decomposition must resolve a spread of up to
its reciprocal; a larger one narrows it, so the recursion runs deeper
and every level adds the rounding of one more reduction. Measured on the
rank identity of
[`additive_check_tol()`](https://statmodels7.github.io/penalties7/reference/additive_check_tol.md),
worst case over nine shapes and 200 parameter vectors each drawn over
the range a fit visits:

|  |  |  |  |  |  |  |
|----|----|----|----|----|----|----|
| `d_tol` | \\\epsilon^{0.25}\\ | \\\epsilon^{0.30}\\ | \\\epsilon^{0.35}\\ | \\\epsilon^{0.40}\\ | \\\epsilon^{0.50}\\ | \\\epsilon^{0.60}\\ |
| worst | 2.8e-01 | 2.8e-01 | 5.3e-05 | 5.3e-05 | 1.4e-03 | 2.2e-03 |

mgcv's `gam.reparam` uses \\\epsilon^{0.3}\\ for the same job. It is
measurably worse on the structures this package builds, which is what
the sweep is for: a constant tuned on another package's penalties is
tuned on another package's spreads. The one-sided normalization an
anisotropic `te()` applies to its margins leaves them incommensurable by
as much as `7.5e+04` before any parameter is estimated, so the spreads
reached here are not the spreads reached there.

## See also

[`additive_stable()`](https://statmodels7.github.io/penalties7/reference/additive_stable.md),
which reads it.
