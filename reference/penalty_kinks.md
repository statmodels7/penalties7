# Points Where a Penalty Is Not Smooth

Returns the values of \\t = D\beta\\ at which some derivative of
\\\rho\\ is discontinuous, so that a numerical reference straddling one
of them would measure the break instead of the formula.
[`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
places its grids clear of these points, and a fitting routine can read
them to determine whether a block can be passed to a gradient method.

## Usage

``` r
penalty_kinks(pen, theta, ...)
```

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object of any branch.

- theta:

  A named list of hyperparameter values, or a named numeric vector
  carrying the same, holding every name in `pen@params`. The argument is
  required even for a smooth penalty, whose result does not depend on
  it.

- ...:

  Passed to methods. The shipped methods do not read it.

## Value

A numeric vector of the points, unsorted and possibly empty, in the
units of \\t = D\beta\\. The result is `numeric(0)` for a penalty that
is smooth everywhere.

## Kinks and the other breaks

A kink in the strict sense is a point where the value is continuous and
the first derivative jumps: \\t = 0\\ for the lasso, the elastic net,
SCAD and MCP, where the subdifferential opens into an interval and a
coefficient can be held exactly at zero. The function returns the larger
set of points where any derivative breaks, because a method that uses
the second derivative needs those points as well. At \\\lambda = 1\\,
\\a = 3.7\\ and \\\gamma = 3\\, with unit curvature, the points are:

|  |  |  |
|----|----|----|
| penalty | returned | what breaks there |
| quadratic, ridge, structured, additive | `numeric(0)` | nothing, the value is smooth |
| separable with a smooth parent (Gaussian, Student t) | `numeric(0)` | nothing |
| lasso, elastic net | `0` | the first derivative |
| SCAD | `0, -1, 1, -3.7, 3.7` | the first at 0, the second at \\\pm\lambda\\ and \\\pm a\lambda\\ |
| MCP | `0, -3, 3` | the first at 0, the second at \\\pm\gamma\lambda\\ |

SCAD's first derivative is continuous at \\\pm\lambda\\, where both
pieces give \\\lambda\\, while its second derivative jumps from \\0\\ to
\\-1/(a-1)\\; at \\\pm a\lambda\\ the second derivative returns to
\\0\\. MCP's first derivative is continuous at \\\pm\gamma\lambda\\, and
its second derivative jumps there from \\-1/\gamma\\ to \\0\\.

## Dependence on the hyperparameters

For a separable penalty the points are fixed when the penalty is built,
from the parameters of the parent that
[`distributions7::fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.html)
holds at a value: the lasso has its kink at \\0\\, and a Laplace parent
whose location is fixed at \\2\\ has its kink at \\2\\. For SCAD and MCP
only \\0\\ is fixed. Their other points are at \\\pm\lambda/c_j\\ and
\\\pm a\lambda/c_j\\ for SCAD and at \\\pm\gamma\lambda/c_j\\ for MCP,
where \\c_j\\ is the curvature of coordinate \\j\\ (one by default), so
the result depends on `theta` and the argument is required.

## See also

[`penalty_gradient()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
for what the first derivative returns at a kink,
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
for the operator that steps over one,
[`has_prox()`](https://statmodels7.github.io/penalties7/reference/has_prox.md)
for whether the penalty has such an operator,
[`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
for the caller that uses this to place a grid.

## Examples

``` r
# Smooth penalties return nothing.
penalty_kinks(quadratic_penalty(diag(2)), list(lambda = 1))
#> numeric(0)

# The lasso breaks only at the origin.
penalty_kinks(lasso_penalty(), list(lambda = 1))
#> [1] 0

# SCAD returns five points: the kink at zero and the two pairs where its
# second derivative changes branch.
penalty_kinks(scad_penalty(), list(lambda = 1, a = 3.7))
#> [1]  0.0 -1.0  1.0 -3.7  3.7

# The points other than zero move with the hyperparameters.
penalty_kinks(scad_penalty(), list(lambda = 2, a = 3.7))
#> [1]  0.0 -2.0  2.0 -7.4  7.4
```
