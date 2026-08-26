# The Proximal Operator as a Piecewise Linear Table

Returns the scalar proximal operator of a separable penalty as an odd
piecewise linear map: three matrices of cuts, slopes and intercepts from
which a compiled loop can apply the operator without knowing which
family it came from. Returns `NULL` for a penalty whose operator has no
such description, and a caller that gets `NULL` uses
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
itself.

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

  A numeric vector of step lengths, one per coefficient, recycled from a
  single value. Every entry must be positive, and **that is not
  checked**: a negative step produces a negative cut and a table nothing
  will reject.
  [`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md),
  which takes a single step, does check it.

- ...:

  Passed to methods. No shipped method reads it.

## Value

A list of three matrices, `cut`, `slope` and `icept`, each with
`pen@n_coef` rows and one column per piece, the last column of `cut`
being `Inf`. `NULL` where the operator has no piecewise linear
description.

## What the table says

The operator of a separable penalty acts one coordinate at a time and is
odd, so it is determined by what it does to \\\lvert u\rvert\\: on the
\\k\\-th interval, \\\lvert u\rvert \le\\ `cut[j, k]`,

\$\$\mathrm{prox}(u) = \mathrm{sign}(u)\\(a\_{jk}\lvert u\rvert +
b\_{jk}).\$\$

Every closed form the package carries has that shape. At \\\lambda =
1.5\\ for the first three, \\\lambda = 1\\, \\a = 3.7\\, \\\gamma = 3\\
for the last two, and a step of `0.5`:

|                        |        |                            |
|------------------------|--------|----------------------------|
| penalty                | pieces | cuts                       |
| Gaussian prior (ridge) | 1      | `Inf`                      |
| Laplace prior (lasso)  | 2      | `0.75`, `Inf`              |
| elastic net            | 2      | `0.45`, `Inf`              |
| MCP                    | 3      | `0.5`, `3`, `Inf`          |
| SCAD                   | 4      | `0.5`, `1.5`, `3.7`, `Inf` |

The last cut is always `Inf`. Applied with
[`prox_apply()`](https://statmodels7.github.io/penalties7/reference/prox_apply.md)
and compared against
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
at every cut, on both sides of every cut, and over a grid across the
whole range, all five agree to `8.9e-16` or better.

## Why a table rather than the operator

A coordinate descent applies the operator once per coordinate per sweep,
at a point that moves every time, so a compiled loop calling back into R
for it would spend its gain on the calls. Passing the numbers lets the
loop stay compiled while the penalty keeps the mathematics: the kernel
evaluates any map of this shape and names no family.

## Why the step is a vector

In a coordinate descent the step of coordinate \\j\\ is \\1/v_j\\ with
\\v_j = \sum_i w_i x\_{ij}^2\\, which does not move while the working
weights are held. The whole table is therefore built once per weighted
least squares iteration and every sweep reads it. Note the asymmetry
with
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md),
which takes a **single** step.

## Under a diagonal map

Standardization is a diagonal \\D\\, under which a separable penalty
stays separable and the table survives. The builder is called at the
step \\t_j d_j^2\\, and the resulting cuts and intercepts are divided by
\\\lvert d_j\rvert\\ while the slopes do not move, the slope multiplying
a point that was scaled and then divided back. The operator is odd, so
only the magnitude of \\d\\ enters. For a lasso this leaves the cut at
\\t_j d_j \lambda\\: at \\d = (0.5, 2, 3)\\, \\t = 0.4\\ and \\\lambda =
1.2\\ the cuts are `0.24`, `0.96`, `1.44`.

The convexity condition of SCAD and MCP is tested on the scaled step, so
it tightens to \\t \< (a-1)/d_j^2\\ and \\t \< \gamma/d_j^2\\: a
standardized penalty takes shorter steps.

## What has no table

`NULL` comes back for a quadratic or structured penalty, whose operator
is one linear solve over every coordinate at once; for a separable
penalty under a map that is not diagonal, which is the generalized-lasso
problem; for one whose parent is read blockwise, whose coordinates do
not separate; for one whose operator is a root rather than a formula,
such as the heavy-tailed Student t prior; for one whose parent is not
centered where the quadratic pull is; and for SCAD or MCP at a step
where the operator is set-valued, that is `step >= a - 1` or
`step >= gamma`.

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

# The table reproduces the operator exactly. The first point sits on its
# own cut and is set to zero; the second survives, shrunk by t * lambda.
u <- c(0.75, 4)
prox_apply(sp, u)
#> [1] 0 1
c(penalty_prox(lasso_penalty(n_coef = 1), u[1], 0.5, list(lambda = 1.5)),
  penalty_prox(lasso_penalty(n_coef = 1), u[2], 2, list(lambda = 1.5)))
#> [1] 0 1

# SCAD needs four pieces, and has none where its operator is set-valued.
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
