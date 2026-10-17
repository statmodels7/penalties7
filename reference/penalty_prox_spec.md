# The Proximal Operator as a Piecewise Linear Table

Returns the scalar proximal operator of a separable penalty as an odd
piecewise linear map: three matrices of cuts, slopes and intercepts from
which a compiled loop can apply the operator without reference to the
family. Returns `NULL` for a penalty whose operator has no such
description; where such a penalty has an operator, it is applied with
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md).

## Usage

``` r
penalty_prox_spec(pen, theta, step, ...)
```

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object of any branch.

- theta:

  A named list of hyperparameter values, or a named numeric vector
  carrying the same, holding every name in `pen@params`.

- step:

  A numeric vector of step lengths, either a single value used for every
  coefficient or one value per coefficient. Every entry must be finite
  and positive, and any other length or value is rejected with an error
  naming `step`.

- ...:

  Passed to methods. The shipped methods do not read it.

## Value

A list of three matrices, `cut`, `slope` and `icept`, each with
`pen@n_coef` rows and one column per piece. The last column of `cut` is
`Inf`. The value is `NULL` where the operator has no piecewise linear
description.

## The piecewise linear map

The operator of a separable penalty acts one coordinate at a time and is
odd, so it is determined by its action on \\\lvert u\rvert\\: on the
\\k\\-th interval, the first for which \\\lvert u\rvert \le\\
`cut[j, k]`,

\$\$\mathrm{prox}(u) = \mathrm{sign}(u)\\(a\_{jk}\lvert u\rvert +
b\_{jk}),\$\$

with \\a\_{jk}\\ = `slope[j, k]` and \\b\_{jk}\\ = `icept[j, k]`. Every
closed-form operator of a separable penalty in the package has this
shape. The table below is for a step of `0.5`, with \\\sigma = 0.8\\ for
the Gaussian prior (whose single cut does not depend on \\\sigma\\),
\\\lambda = 1.5\\ for the Laplace prior, \\\lambda = 1.5\\ and \\\alpha
= 0.6\\ for the elastic net, \\\lambda = 1\\ and \\\gamma = 3\\ for MCP,
and \\\lambda = 1\\ and \\a = 3.7\\ for SCAD:

|  |  |  |
|----|----|----|
| penalty | pieces | cuts |
| Gaussian prior (`gaussian1` parent) | 1 | `Inf` |
| Laplace prior ([`lasso_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md)) | 2 | `0.75`, `Inf` |
| elastic net ([`elasticnet_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md)) | 2 | `0.45`, `Inf` |
| MCP | 3 | `0.5`, `3`, `Inf` |
| SCAD | 4 | `0.5`, `1.5`, `3.7`, `Inf` |

The last cut is always `Inf`.
[`prox_apply()`](https://statmodels7.github.io/penalties7/reference/prox_apply.md)
evaluates a table, and the tests of the package compare it with
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
at every cut, on both sides of every cut and over a grid of points.

## Accuracy near the step limit of SCAD and MCP

The table of SCAD and MCP is ill-conditioned as the step approaches its
limit, \\a - 1\\ for SCAD and \\\gamma\\ for MCP. On the tapering piece
the slope is \\1/(1 - t/(a-1))\\ for SCAD and \\1/(1 - t/\gamma)\\ for
MCP, which grows without bound as the step approaches the limit, and the
intercept grows with it and has the opposite sign, so \\a\_{jk}\lvert
u\rvert + b\_{jk}\\ is the difference of two large numbers that nearly
cancel. The table is therefore less accurate there than
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md),
which evaluates the same piece in a form without this cancellation. With
a curvature or a diagonal map the limit applies to the scaled step \\t_j
c_j d_j^2\\ described below.

## Why a table rather than the operator

A coordinate descent applies the operator once per coordinate per sweep,
at a point that moves every time, so a compiled loop that called back
into R for each application would lose the speed gained by compiling it.
Passing the numbers keeps the loop compiled. The table holds the
mathematics of the penalty, and the compiled kernel evaluates any map of
this shape without reference to the family.

## Why the step is a vector

In a coordinate descent the step of coordinate \\j\\ is \\1/v_j\\ with
\\v_j = \sum_i w_i x\_{ij}^2\\, which does not change while the working
weights are held fixed. The whole table is therefore built once per
weighted least squares iteration and read by every sweep.
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md),
by contrast, takes a single step.

## Under a diagonal map

Standardization is a diagonal \\D\\, under which a separable penalty
stays separable and the table is still available. The builder is called
at the step \\t_j d_j^2\\, and the cuts and intercepts that it returns
are divided by \\\lvert d_j\rvert\\, while the slopes are unchanged
because each slope multiplies a point that was scaled and then divided
back. The operator is odd, so only the magnitude of \\d\\ enters. For a
lasso the cut is \\t_j \lvert d_j\rvert \lambda\\ (with \\d = (0.5, 2,
3)\\, \\t = 0.4\\ and \\\lambda = 1.2\\ the cuts are `0.24`, `0.96` and
`1.44`).

The step limit of SCAD and MCP is tested on the scaled step. With
\\c_j\\ the curvature of the coordinate (the `curv` property, one when
it is empty), the table exists only while \\t_j c_j d_j^2 \< a - 1\\ for
SCAD and \\t_j c_j d_j^2 \< \gamma\\ for MCP in every coordinate, so a
standardized penalty takes shorter steps.

A diagonal map is recognized by its class through
[`map_diagonal()`](https://statmodels7.github.io/penalties7/reference/map_diagonal.md).
[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md)
(and therefore
[`lasso_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md)
and
[`elasticnet_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md)),
[`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
and
[`mcp_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
store a square map whose off-diagonal entries are all zero as a
[`Matrix::Diagonal()`](https://rdrr.io/pkg/Matrix/man/Diagonal.html)
object, so a diagonal base matrix gives a table. A diagonal map with a
zero entry gives `NULL`.

## Penalties without a table

`NULL` is returned for a quadratic or structured penalty, whose operator
is one linear solve over every coordinate at once, and for an additive
penalty. Among separable penalties, `NULL` is returned under a map that
is not recognized as diagonal, which is the generalized-lasso problem;
for a parent read blockwise, whose coordinates do not separate; for a
parent other than the Gaussian, Laplace and elastic-net families, whose
operator has no closed form (the Student t prior, for example, has an
operator computed by root finding); for a parent whose gradient at zero
is not zero, that is, a parent not centered at zero; and for SCAD or MCP
when the scaled step \\t_j c_j d_j^2\\ of some coordinate reaches \\a -
1\\ or \\\gamma\\, where the operator is set-valued.

## See also

[`prox_apply()`](https://statmodels7.github.io/penalties7/reference/prox_apply.md)
to evaluate the table,
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
for the operator itself,
[`has_prox()`](https://statmodels7.github.io/penalties7/reference/has_prox.md)
for whether there is one at all,
[`spec_diag()`](https://statmodels7.github.io/penalties7/reference/spec_diag.md)
for the diagonal-map transport.

## Examples

``` r
# A lasso at two different steps: the threshold is t * lambda.
pen <- lasso_penalty(n_coef = 2L)
sp <- penalty_prox_spec(pen, list(lambda = 1.5), step = c(0.5, 2))
sp$cut
#>      [,1] [,2]
#> [1,] 0.75  Inf
#> [2,] 3.00  Inf

# The table reproduces the operator. The first point lies on its cut and
# is set to zero; the second is shrunk by t * lambda.
u <- c(0.75, 4)
prox_apply(sp, u)
#> [1] 0 1
c(penalty_prox(lasso_penalty(n_coef = 1), u[1], 0.5, list(lambda = 1.5)),
  penalty_prox(lasso_penalty(n_coef = 1), u[2], 2, list(lambda = 1.5)))
#> [1] 0 1

# SCAD has four pieces, and no table at a step where its operator is
# set-valued.
penalty_prox_spec(scad_penalty(), list(lambda = 1, a = 3.7),
                  step = 0.5)$cut
#>      [,1] [,2] [,3] [,4]
#> [1,]  0.5  1.5  3.7  Inf
penalty_prox_spec(scad_penalty(), list(lambda = 1, a = 3.7), step = 2.7)
#> NULL

# A quadratic penalty has a solve, not a coordinatewise map.
penalty_prox_spec(quadratic_penalty(diag(2)), list(lambda = 1), step = 0.5)
#> NULL
```
