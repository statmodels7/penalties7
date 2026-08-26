# Construct the SCAD and MCP Penalties

The two non-convex selection penalties that shrink small coefficients
like the lasso and leave large ones alone. Both are defined by their
derivative on \\t \ge 0\\ and extended evenly, so `scad_penalty()` and
`mcp_penalty()` build objects whose value is the closed antiderivative
anchored at \\\rho(0) = 0\\.

## Usage

``` r
scad_penalty(
  map = NULL,
  n_coef = 1L,
  link_lambda = linkfunctions7::log_link(),
  link_a = linkfunctions7::bounded_link(lwr = 2)
)

mcp_penalty(
  map = NULL,
  n_coef = 1L,
  link_lambda = linkfunctions7::log_link(),
  link_gamma = linkfunctions7::bounded_link(lwr = 1)
)
```

## Arguments

- map:

  The matrix \\D\\, with one column per coefficient, or `NULL` (the
  default) for the identity. Given a map, `n_coef` is taken from its
  column count and the argument is ignored.

- n_coef:

  The number of coefficients when `map` is `NULL`. A single whole
  number, `1L` by default.

- link_lambda:

  The linkfunctions7 link carrying \\\lambda\\ onto the whole real line.
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  by default, \\\lambda\\ being positive.

- link_a:

  The link carrying SCAD's shape parameter, bounded below at 2.
  `linkfunctions7::bounded_link(lwr = 2)` by default. `scad_penalty()`
  only.

- link_gamma:

  The link carrying MCP's shape parameter, bounded below at

  1.  `linkfunctions7::bounded_link(lwr = 1)` by default.
      `mcp_penalty()` only.

## Value

`scad_penalty()` a
[`ScadPenalty()`](https://statmodels7.github.io/penalties7/reference/ScadPenalty.md)
object with hyperparameters `lambda` on \\(0, \infty)\\ and `a` on \\(2,
\infty)\\. `mcp_penalty()` an
[`McpPenalty()`](https://statmodels7.github.io/penalties7/reference/ScadPenalty.md)
object with `lambda` on \\(0, \infty)\\ and `gamma` on \\(1, \infty)\\.

## The two derivatives

SCAD, after Fan and Li (2001), holds the lasso's slope up to
\\\lambda\\, tapers it linearly to zero over the next \\(a-1)\lambda\\,
and stops:

\$\$\rho'(t) = \lambda \\\\ (t \le \lambda), \qquad \rho'(t) =
\frac{a\lambda - t}{a - 1} \\\\ (\lambda \< t \le a\lambda), \qquad
\rho'(t) = 0 \\\\ (t \> a\lambda).\$\$

MCP, after Zhang (2010), starts tapering at once:

\$\$\rho'(t) = \left(\lambda - t/\gamma\right)\_{+}.\$\$

Integrating gives the values, which saturate at \\(a+1)\lambda^2/2\\ and
\\\gamma\lambda^2/2\\. Both agree with a quadrature of the published
derivative to `4e-15` or better at every argument, in each of the three
regions and at the boundaries between them.

## What the shape parameter buys

Near zero both behave like a lasso at rate \\\lambda\\, so they
threshold and produce exact zeros. Far from zero both are flat, so a
large coefficient is estimated without shrinkage, so the lasso's bias is
gone. The shape parameter says how quickly the transition happens: \\a
\to \infty\\ and \\\gamma \to \infty\\ recover the lasso, while small
values approach hard thresholding. The literature's defaults are \\a =
3.7\\ and \\\gamma = 3\\, and neither constructor supplies one, the
value being the caller's to set or a path's to sweep.

## Where they are not smooth

Both have a kink at zero, and both change branch again where the taper
begins or ends.
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
returns all of them: `0`, \\\pm\lambda\\ and \\\pm a\lambda\\ for SCAD,
`0` and \\\pm\gamma\lambda\\ for MCP. At the outer points the first
derivative is continuous and the second jumps.
[`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
asks the object and keeps its grids away from all of them.

## What they are not

Both are improper: \\\rho\\ is bounded, so \\\exp(-\rho)\\ is not a
density at any constant.
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
is `FALSE`,
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
is `FALSE`, and
[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
and its three siblings reject. A marginal criterion cannot reach these
hyperparameters; a path or cross-validation can.

## References

Fan, J. and Li, R. (2001). Variable selection via nonconcave penalized
likelihood and its oracle properties. *Journal of the American
Statistical Association* **96**, 1348-1360.

Zhang, C.-H. (2010). Nearly unbiased variable selection under minimax
concave penalty. *Annals of Statistics* **38**, 894-942.

## See also

[`lasso_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md)
for the convex penalty these improve on,
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
for the closed operator each has over its convex region,
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
for where they change branch,
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
for why no marginal criterion reaches them.

## Examples

``` r
# SCAD at the literature's shape. Small coefficients pay the lasso's rate,
# large ones pay a constant.
pen <- scad_penalty(n_coef = 3)
th <- list(lambda = 1, a = 3.7)
penalty_value(pen, c(0.5, 2, 5), th)
#> [1] 4.664815

# The value saturates at (a + 1) lambda^2 / 2 per coefficient.
penalty_value(scad_penalty(n_coef = 1), 100, th)
#> [1] 2.35
(3.7 + 1) * 1^2 / 2
#> [1] 2.35

# So the gradient of a large coefficient is exactly zero: no shrinkage.
penalty_gradient(pen, c(0.5, 2, 5), th)
#> [1] 1.0000000 0.6296296 0.0000000

# MCP saturates sooner, at gamma lambda^2 / 2.
mcp <- mcp_penalty(n_coef = 3)
penalty_value(mcp, c(0.5, 2, 5), list(lambda = 1, gamma = 3))
#> [1] 3.291667
penalty_gradient(mcp, c(0.5, 2, 5), list(lambda = 1, gamma = 3))
#> [1] 0.8333333 0.3333333 0.0000000

# Both are bounded, so neither is a density.
c(scad = is_proper(pen), mcp = is_proper(mcp))
#>  scad   mcp 
#> FALSE FALSE 

# Letting the shape run recovers the lasso's linear penalty.
sapply(c(3.7, 50, 1000),
       function(a) penalty_value(scad_penalty(n_coef = 1), 2,
                                 list(lambda = 1, a = a)))
#> [1] 1.814815 1.989796 1.999499
```
