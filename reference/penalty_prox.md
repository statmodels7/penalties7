# Proximal Operator of a Penalty

Returns the point that minimizes the penalty plus a quadratic pull
towards \\v\\,

\$\$\mathrm{prox}\_{t\rho}(v) = \arg\min\_{\beta} \left\\
\tfrac{1}{2t}\lVert \beta - v \rVert^{2} + \rho(\beta;\theta)
\right\\,\$\$

A proximal gradient method applies the operator after each gradient
step. With it, a penalty with a kink is minimized without
differentiating it at the kink, and a coefficient can be set exactly to
zero, where a gradient method only brings it close to zero.

## Usage

``` r
penalty_prox(penalty, v, step, theta, ...)
```

## Arguments

- penalty:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.
  [`has_prox()`](https://statmodels7.github.io/penalties7/reference/has_prox.md)
  tests whether it has an operator.

- v:

  A numeric vector of length `penalty@n_coef`, the point of the
  quadratic pull. The vector is coerced with
  [`as.numeric()`](https://rdrr.io/r/base/numeric.html).

- step:

  The step length \\t\\, a single positive finite number. A vector is
  rejected;
  [`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
  takes one step per coefficient.

- theta:

  A named list of hyperparameter values, or a named numeric vector,
  holding every name in `penalty@params`, each a single number.

- ...:

  Passed to methods. The shipped methods do not use it.

## Value

A numeric vector of the same length as `v`.

## The operator of each branch

|  |  |
|----|----|
| branch | operator |
| quadratic, structured | one linear solve, \\(I + tS)^{-1}v\\, with \\S\\ the penalty's Hessian |
| Gaussian prior centered at zero | the shrinkage \\v/(1 + t/\sigma^2)\\ |
| Laplace prior | the soft threshold \\\mathrm{sign}(v)(\lvert v\rvert - t\lambda)\_{+}\\ |
| elastic net | the soft threshold at \\t\lambda\alpha\\, followed by division by \\1 + t\lambda(1-\alpha)\\ |
| SCAD, MCP | the piecewise operators of Fan and Li (2001) and Zhang (2010) |
| any other separable parent | the stationary point of smallest objective, coordinate by coordinate |

Each closed form is the minimizer of the defining objective, and the
page of each method gives its formula.

## The stationary condition, where there is no closed form

A separable penalty whose parent has a differentiable log-density, and
which is not one of the closed-form families, is solved coordinate by
coordinate from the stationary condition

\$\$\frac{\beta - v}{t} = \ell^{(y)}(\beta),\$\$

whose right-hand side is
[`distributions7::distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.html).
For a density that is not log-concave, such as the Student t of
[`heavy_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md),
the condition can have more than one root, and the operator is the root
of smallest objective. The segment between \\v\\ and the mode of the
parent is scanned for changes of sign, each change is refined by
[`uniroot()`](https://rdrr.io/r/stats/uniroot.html) to a relative
accuracy of about `eps^0.75`, and the candidate of smallest objective is
returned. For a unimodal parent the global minimizer lies on that
segment, so the value returned is the global minimizer. For a parent
with more than one mode, the value returned is a stationary point that
need not be the global minimizer;
[`.prox_root()`](https://statmodels7.github.io/penalties7/reference/dot-prox_root.md)
gives the details. A parent that declares a kink and is not one of the
closed-form families is rejected, because the stationary condition does
not hold at a kink and the minimizer of such a penalty can lie there.

## Under a map

The quadratic branch takes any map, because the objective stays
quadratic, and a structured penalty has no map. The separable, SCAD and
MCP branches need the identity map or a diagonal one. Under a diagonal
map with entries \\d_j\\, the operator is the operator without the map
evaluated at \\d_j v_j\\ with the step \\t d_j^2\\, and then divided by
\\d_j\\; in this way the standardization of a column enters the
operator. A general \\D\\ mixes the coordinates and turns the problem
into the generalized-lasso problem, which requires an algorithm of its
own, so such a penalty is rejected with a message naming it. A diagonal
map is recognized by its Matrix class, as
[`map_diagonal()`](https://statmodels7.github.io/penalties7/reference/map_diagonal.md)
describes.

## Calls that signal an error

- a penalty with no operator, such as an
  [`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md),
  with a message naming
  [`has_prox()`](https://statmodels7.github.io/penalties7/reference/has_prox.md);

- a separable, SCAD or MCP penalty under a map that is not a Matrix
  diagonal object with nonzero entries;

- a separable penalty whose parent is read in blocks, because the
  coordinates of a block do not separate;

- a Laplace or elastic-net parent that is not centered at zero, because
  the closed form is written for a parent centered at zero;

- a separable parent that declares a kink and has no closed form;

- a stationary condition that changes sign neither on the segment
  between \\v\\ and the mode of the parent nor in a widened bracket;

- SCAD where \\t c_j \ge a - 1\\ and MCP where \\t c_j \ge \gamma\\ in
  some coordinate, where \\c_j\\ is the `curvature` property of the
  penalty and \\t\\ is multiplied by \\d_j^2\\ under a diagonal map; the
  subproblem is then not convex and the operator is set-valued;

- a `step` that is not a single positive finite number, or a `v` whose
  length is not `penalty@n_coef`.

## Methods

The method on the base class
[`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
signals an error naming the penalty and pointing to
[`has_prox()`](https://statmodels7.github.io/penalties7/reference/has_prox.md);
[`AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/AdditivePenalty.md)
inherits it. The classes with an operator of their own are:

- one linear solve;

- [`StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/StructuredPenalty.md):
  [`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.QuadraticPenalty.md),
  the same linear solve;

- a closed form for the Gaussian, Laplace and elastic-net parents and
  the stationary point of smallest objective for any other parent;

- [`ScadPenalty()`](https://statmodels7.github.io/penalties7/reference/ScadPenalty.md):
  [`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.ScadPenalty.md),
  the piecewise operator of Fan and Li (2001);

- [`McpPenalty()`](https://statmodels7.github.io/penalties7/reference/ScadPenalty.md):
  [`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.McpPenalty.md),
  the piecewise operator of Zhang (2010).

## References

Fan, J. and Li, R. (2001). Variable selection via nonconcave penalized
likelihood and its oracle properties. *Journal of the American
Statistical Association* **96**, 1348-1360.

Zhang, C.-H. (2010). Nearly unbiased variable selection under minimax
concave penalty. *Annals of Statistics* **38**, 894-942.

## See also

[`has_prox()`](https://statmodels7.github.io/penalties7/reference/has_prox.md)
to test whether a penalty has an operator,
[`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
for the same operator as a table that a compiled loop can read,
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
for the points where a gradient method fails,
[`penalty_gradient()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
for the gradient where the penalty is smooth.

## Examples

``` r
v <- c(2, 0.3, -1.4)

# The lasso: a soft threshold, so the middle coordinate becomes exactly 0.
penalty_prox(lasso_penalty(n_coef = 3), v, 1, list(lambda = 1))
#> [1]  1.0  0.0 -0.4

# A ridge: a shrinkage towards zero, which sets no coordinate to zero.
penalty_prox(ridge_penalty(n_coef = 3), v, 1, list(lambda = 1.5))
#> [1]  0.80  0.12 -0.56
v / (1 + 1 * 1.5)
#> [1]  0.80  0.12 -0.56

# The elastic net applies the soft threshold and then the shrinkage.
penalty_prox(elasticnet_penalty(n_coef = 3), v, 1,
             list(lambda = 1, alpha = 0.6))
#> [1]  1.0000000  0.0000000 -0.5714286

# A heavy-tailed prior has no closed form, so the operator is found from
# the stationary condition. With sigma = 1, nu = 4 and step 1 it maps 2 to
# 1, as the ridge with lambda = 1 does, and it moves the coordinate 10 much
# less than that ridge does.
v2 <- c(2, 0.3, 10)
penalty_prox(heavy_penalty(n_coef = 3), v2, 1, list(sigma = 1, nu = 4))
#> [1] 1.0000000 0.1336635 9.4958181
penalty_prox(ridge_penalty(n_coef = 3), v2, 1, list(lambda = 1))
#> [1] 1.00 0.15 5.00

# Each operator minimizes the defining objective. The next call finds the
# minimizer for the first coordinate of the lasso.
rho <- lasso_penalty(n_coef = 3)
obj <- function(b1) {
  x <- v; x[1] <- b1
  0.5 * sum((x - v)^2) / 1 + penalty_value(rho, x, list(lambda = 1))
}
stats::optimize(obj, c(-5, 5), tol = 1e-12)$minimum
#> [1] 1

# An additive penalty has no operator, and the call signals an error.
try(penalty_prox(additive_penalty(list(diag(3), diag(c(1, 1, 0)))),
                 v, 1, list(lambda1 = 1, lambda2 = 1)))
#> Error : 'additive [2 components]' has no proximal operator.
#>   A penalty has one when it is a quadratic or structured penalty, a SCAD
#>   or MCP penalty, or a separable penalty under the identity map or a
#>   diagonal map. has_prox() tests this before a call.
```
