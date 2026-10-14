# Construct the SCAD and MCP Penalties

The two non-convex selection penalties that shrink small coefficients
like the lasso and leave large ones unshrunk. Both are defined by their
derivative on \\t \ge 0\\ and extended evenly, and `scad_penalty()` and
`mcp_penalty()` build objects whose value is the closed antiderivative
anchored at \\\rho(0) = 0\\.

## Usage

``` r
scad_penalty(
  map = NULL,
  n_coef = 1L,
  link_lambda = linkfunctions7::log_link(),
  link_a = linkfunctions7::bounded_link(lwr = 2),
  curv = NULL
)

mcp_penalty(
  map = NULL,
  n_coef = 1L,
  link_lambda = linkfunctions7::log_link(),
  link_gamma = linkfunctions7::bounded_link(lwr = 1),
  curv = NULL
)
```

## Arguments

- map:

  The matrix \\D\\, with one column per coefficient, or `NULL` (the
  default) for the identity. Given a map, `n_coef` is taken from its
  column count and the argument is ignored. A Matrix object keeps its
  storage, and a square map whose off-diagonal entries are all zero is
  stored as a
  [`Matrix::Diagonal()`](https://rdrr.io/pkg/Matrix/man/Diagonal.html),
  so that a diagonal map has a proximal operator whatever its class.

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

- curv:

  `NULL` (the default) for unit curvature, or positive finite numbers,
  one per coordinate of \\D\beta\\ or a single value repeated.

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

SCAD, after Fan and Li (2001), keeps the slope \\\lambda\\ of the lasso
up to \\\lambda\\, tapers it linearly to zero over the next
\\(a-1)\lambda\\, and is zero beyond:

\$\$\rho'(t) = \lambda \\\\ (t \le \lambda), \qquad \rho'(t) =
\frac{a\lambda - t}{a - 1} \\\\ (\lambda \< t \le a\lambda), \qquad
\rho'(t) = 0 \\\\ (t \> a\lambda).\$\$

MCP, after Zhang (2010), tapers from the origin:

\$\$\rho'(t) = \left(\lambda - t/\gamma\right)\_{+}.\$\$

Integrating gives the values, which saturate at \\(a+1)\lambda^2/2\\ for
SCAD and at \\\gamma\lambda^2/2\\ for MCP. These are the forms for unit
curvature; the section on the curvature of the loss gives the form that
`curv` produces.

## The shape parameter

Near zero both penalties behave like a lasso at rate \\\lambda\\, so
they threshold and produce exact zeros. Far from zero both are flat, so
a large coefficient is estimated without shrinkage and without the bias
of the lasso. The shape parameter sets how quickly the transition
happens: \\a \to \infty\\ and \\\gamma \to \infty\\ recover the lasso,
smaller values make the transition more abrupt, and MCP approaches hard
thresholding as \\\gamma \to 1\\. The values used in the literature are
\\a = 3.7\\ and \\\gamma = 3\\. Neither constructor supplies a value,
which the caller sets or a path sweeps.

## Points of non-smoothness

Both penalties have a kink at zero, and both change branch again where
the taper begins or ends.
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
returns all of these points: at unit curvature they are `0`,
\\\pm\lambda\\ and \\\pm a\lambda\\ for SCAD and `0` and
\\\pm\gamma\lambda\\ for MCP, and with `curv` the rate \\\lambda\\ in
them becomes \\\lambda/c_j\\. At the points other than the origin the
first derivative is continuous and the second jumps. When `beta` is not
supplied,
[`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
draws its coefficients away from the points that
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
returns.

## The curvature of the loss

Both penalties are defined on the canonical problem
\\\tfrac12(z-\theta)^2 + \rho(\theta)\\, whose loss has unit curvature
(Fan and Lv 2010), and the value \\a = 3.7\\ of Fan and Li (2001) is
suggested for that problem. A likelihood has curvature \\c_j\\ in
coordinate \\j\\, and with `curv` the penalty is \\\sum_j
c_j\\\rho(t_j;\\ \lambda/c_j)\\. The slope at zero stays \\\lambda\\,
the knee moves to \\a\lambda/c_j\\ (\\\gamma\lambda/c_j\\ for MCP), and
the ratio of the knee to \\\lambda/c_j\\ is the shape parameter in every
coordinate. Coordinate \\j\\ saturates at \\(a+1)\lambda^2/(2c_j)\\ for
SCAD and \\\gamma\lambda^2/(2c_j)\\ for MCP, and its second derivative
over the taper is \\-c_j/(a-1)\\ for SCAD and \\-c_j/\gamma\\ for MCP.
The proximal operator with step \\s\\ is the unscaled operator at step
\\s c_j\\ and rate \\\lambda/c_j\\, so it is single-valued while \\s c_j
\< a - 1\\ for SCAD and \\s c_j \< \gamma\\ for MCP. At the step \\s =
1/c_j\\ these conditions are \\a \> 2\\ and \\\gamma \> 1\\ for any
data. By default `curv` is `NULL`, no curvature is stored, and every
coordinate has \\c_j = 1\\.

## Improper penalties

Both penalties are improper: \\\rho\\ is bounded, so \\\exp(-\rho)\\ is
not a density at any constant.
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
and
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
return `FALSE`, and
[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
and the three quantities beside it signal an error. A marginal criterion
therefore cannot be used to choose these hyperparameters, which are
chosen along a path or by cross-validation.

## References

Fan, J. and Li, R. (2001). Variable selection via nonconcave penalized
likelihood and its oracle properties. *Journal of the American
Statistical Association* **96**, 1348-1360.

Fan, J. and Lv, J. (2010). A selective overview of variable selection in
high dimensional feature space. *Statistica Sinica* **20**, 101-148.

Zhang, C.-H. (2010). Nearly unbiased variable selection under minimax
concave penalty. *Annals of Statistics* **38**, 894-942.

## See also

[`lasso_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md)
for the convex penalty that both modify,
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
for the closed proximal operator, defined for steps below \\a - 1\\ for
SCAD and below \\\gamma\\ for MCP,
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
for the points where they change branch, and
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
for the property that excludes a marginal criterion.

## Examples

``` r
# SCAD at the shape used in the literature. A small coefficient is
# penalized at the rate of the lasso, a large one by a constant.
pen <- scad_penalty(n_coef = 3)
th <- list(lambda = 1, a = 3.7)
penalty_value(pen, c(0.5, 2, 5), th)
#> [1] 4.664815

# The value saturates at (a + 1) lambda^2 / 2 per coefficient.
penalty_value(scad_penalty(n_coef = 1), 100, th)
#> [1] 2.35
(3.7 + 1) * 1^2 / 2
#> [1] 2.35

# The gradient at a large coefficient is exactly zero, so that coefficient
# is not shrunk.
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

# As the shape grows, the value at 2 approaches the lasso's value 2.
sapply(c(3.7, 50, 1000),
       function(a) penalty_value(scad_penalty(n_coef = 1), 2,
                                 list(lambda = 1, a = a)))
#> [1] 1.814815 1.989796 1.999499

# With curv = 2 the slope at zero stays lambda and the value saturates at
# (a + 1) lambda^2 / (2 * 2).
penalty_value(scad_penalty(n_coef = 1, curv = 2), 100, th)
#> [1] 1.175
(3.7 + 1) * 1^2 / (2 * 2)
#> [1] 1.175
```
