# The Proximal Operator as a Piecewise Linear Table

Describes the scalar proximal operator of a separable penalty as an odd
piecewise linear map, so that a compiled loop can apply it without
knowing which family it came from.

## Usage

``` r
penalty_prox_spec(pen, theta, step, ...)
```

## Arguments

- pen:

  A
  [`penalty`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- theta:

  A named list of hyperparameter values.

- step:

  A numeric vector of step lengths, one per coefficient.

- ...:

  Passed to methods.

## Value

A list with `cut`, `slope` and `icept`, each a matrix with one row per
coefficient and one column per piece, the last cut being `Inf`; or
`NULL` where the operator has no such description.

## Details

**What the table says.** The operator of a separable penalty acts one
coordinate at a time and is odd, so it is determined by what it does to
\\\|u\|\\: on the \\k\\-th interval, \\\|u\| \le\\ `cut[j, k]`,
\$\$\mathrm{prox}(u) = \mathrm{sign}(u)\\(a\_{jk}\|u\| + b\_{jk}).\$\$
Every closed form the package carries has this shape: the soft threshold
is two pieces, the elastic net two, MCP three and SCAD four, and a
Gaussian prior is the single piece \\u/(1 + t/\sigma^2)\\.

**Why a table rather than the operator.** A coordinate descent applies
the operator once per coordinate per sweep, at a point that moves every
time, so a compiled loop calling back into R for it would spend its gain
on the calls – the property that decided how far the score-driven filter
of modelterms7 could be ported. Passing the numbers instead lets the
loop stay compiled while the penalty keeps the mathematics: the kernel
evaluates any map of this shape and names no family.

**Why the step is a vector.** In a coordinate descent the step of
coordinate \\j\\ is \\1/v_j\\ with \\v_j = \sum_i w_i x\_{ij}^2\\, which
does not move while the working weights are held, so the whole table is
built once per weighted least squares iteration and the sweeps read it.

**A diagonal map.** Standardization is a diagonal \\D\\, under which a
separable penalty stays separable and the table survives: reading it at
\\d_j v_j\\ with the step \\t_j d_j^2\\ and dividing back gives the cuts
and the intercepts divided by \\\|d_j\|\\ and the slopes unchanged. The
convexity condition of SCAD and MCP is tested on the scaled step, so it
becomes \\t \< (a-1)/d_j^2\\ and \\t \< \gamma/d_j^2\\.

**What has no table.** A quadratic penalty under a general matrix is not
separable and returns `NULL`, as does a separable penalty under a map
that is not diagonal, one whose operator is a root rather than a
formula, and one whose parent is not centered where the quadratic pull
is. A caller that gets `NULL` uses
[`penalty_prox`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
itself.

## See also

[`penalty_prox`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md),
[`has_prox`](https://statmodels7.github.io/penalties7/reference/has_prox.md)

## Examples

``` r
pen <- lasso_penalty(n_coef = 2L)
penalty_prox_spec(pen, list(lambda = 1.5), step = c(0.5, 2))
#> $cut
#>      [,1] [,2]
#> [1,] 0.75  Inf
#> [2,] 3.00  Inf
#> 
#> $slope
#>      [,1] [,2]
#> [1,]    0    1
#> [2,]    0    1
#> 
#> $icept
#>      [,1]  [,2]
#> [1,]    0 -0.75
#> [2,]    0 -3.00
#> 
```
