# Where a Penalty Stops Being Smooth

Returns the values of \\t = D\beta\\ at which some derivative of
\\\rho\\ is discontinuous, so that a numerical reference straddling one
of them measures the break instead of the formula.
[`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
places its grids clear of these points, and a solver consults them to
decide whether a block can go to a gradient method at all.

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
  carrying the same, holding every name in `pen@params`. Required even
  for a smooth penalty, whose answer does not depend on it.

- ...:

  Passed to methods. No shipped method reads it.

## Value

A numeric vector of the points, unsorted and possibly empty, in the
units of \\t = D\beta\\. `numeric(0)` for a penalty that is smooth
everywhere.

## The set is wider than the kinks alone

A kink proper is a point where the value is continuous and the first
derivative jumps: \\t = 0\\ for the lasso, the elastic net, SCAD and
MCP, where the subdifferential opens into an interval and a coefficient
can be held exactly at zero. What is returned is the larger set of
points where any derivative breaks, because that is what a consumer of
the second derivative needs as well. Measured at \\\lambda = 1, a = 3.7,
\gamma = 3\\:

|  |  |  |
|----|----|----|
| penalty | returned | what breaks there |
| quadratic, ridge, structured, additive | `numeric(0)` | nothing, the value is smooth |
| lasso, elastic net | `0` | the first derivative |
| SCAD | `0, -1, 1, -3.7, 3.7` | the first at 0, the second at \\\pm\lambda\\ and \\\pm a\lambda\\ |
| MCP | `0, -3, 3` | the first at 0, the second at \\\pm\gamma\lambda\\ |

At \\\pm\lambda\\ SCAD's derivative is continuous, both branches giving
\\\lambda\\, while its second derivative jumps from \\0\\ to
\\-1/(a-1)\\. The same holds for MCP at \\\pm\gamma\lambda\\.

## The set moves with the hyperparameters

Only \\0\\ is fixed. SCAD's outer points are at \\\pm\lambda\\ and \\\pm
a\lambda\\ and MCP's at \\\pm\gamma\lambda\\, so the answer depends on
`theta` and the argument is not optional. This is also what a
hyperparameter path is walked over: the size of the kink, rather than
the hyperparameter itself, is the quantity that scales the same way
across branches.

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

# The lasso breaks at the origin alone.
penalty_kinks(lasso_penalty(), list(lambda = 1))
#> [1] 0

# SCAD returns five points: the kink at zero and the two pairs where its
# second derivative changes branch.
penalty_kinks(scad_penalty(), list(lambda = 1, a = 3.7))
#> [1]  0.0 -1.0  1.0 -3.7  3.7

# Which move with the hyperparameters.
penalty_kinks(scad_penalty(), list(lambda = 2, a = 3.7))
#> [1]  0.0 -2.0  2.0 -7.4  7.4
```
