# Proximal Operator of a Penalty

Returns the point that minimizes the penalty plus a quadratic pull
towards \\v\\,

\$\$\mathrm{prox}\_{t\rho}(v) = \arg\min\_{\beta} \left\\
\tfrac{1}{2t}\lVert \beta - v \rVert^{2} + \rho(\beta;\theta)
\right\\,\$\$

A proximal gradient method applies it after each gradient step. It is
the operation that lets a penalty with a kink be minimized without
differencing it, and the one that sets a coefficient exactly to zero
where a gradient method only makes it small.

## Usage

``` r
penalty_prox(pen, v, step, theta, ...)
```

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object. Ask
  [`has_prox()`](https://statmodels7.github.io/penalties7/reference/has_prox.md)
  first if you do not know the branch.

- v:

  A numeric vector of length `pen@n_coef`, the point pulled towards.
  Coerced with [`as.numeric()`](https://rdrr.io/r/base/numeric.html).

- step:

  The step length \\t\\, a single positive number. A vector is rejected;
  [`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
  is the interface that takes one per coefficient.

- theta:

  A named list of hyperparameter values, or a named numeric vector
  carrying the same, holding every name in `pen@params`.

- ...:

  Passed to methods. No shipped method reads it.

## Value

A numeric vector of the same length as `v`.

## The four routes

|  |  |
|----|----|
| branch | operator |
| quadratic, structured | one linear solve, \\(I + tS)^{-1}v\\ |
| Gaussian prior | the shrinkage \\v/(1 + t/\sigma^2)\\ |
| Laplace prior | the soft threshold \\\mathrm{sign}(v)(\lvert v\rvert - t\lambda)\_{+}\\ |
| elastic net | the soft threshold followed by the shrinkage |
| SCAD, MCP | the closed piecewise operators of their papers |
| any other separable parent | a coordinatewise root, found by bisection |

Checked against a direct minimization of the defining objective,
coordinate by coordinate on a three-vector at a step of 1, all eight
shipped penalties agree to `2.7e-08` or better, which is the reference
minimizer's own accuracy.

## The root, where there is no formula

Any separable penalty whose parent has a differentiable log-density is
solved from the stationary condition

\$\$\frac{\beta - v}{t} = \ell^{(y)}(\beta),\$\$

with the right-hand side
[`distributions7::distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.html),
closed form for every continuous family. Log-concavity of the density
makes the left side minus the right side strictly increasing, so the
root is unique; the bracket is widened up to sixty times and then
[`uniroot()`](https://rdrr.io/r/stats/uniroot.html) finds it to
`eps^0.75`. A parent that declares a kink and is not one of the
closed-form families is rejected instead, a smooth root being the wrong
instrument there.

## Under a map

The quadratic and structured branches take **any** map, the objective
staying quadratic. The separable branches need the identity map or a
diagonal one, under which the operator is the unmapped one read at
\\dv\\ with the step \\td^2\\, divided back by \\d\\: that is what
standardization comes to. A general \\D\\ mixes coordinates and makes
this the generalized-lasso problem, which needs an algorithm of its own,
so it is rejected by name.

## What rejects

- a penalty with no operator at all,
  [`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
  being the one that ships, with a message naming
  [`has_prox()`](https://statmodels7.github.io/penalties7/reference/has_prox.md);

- a separable penalty under a map that is not diagonal;

- a separable penalty whose parent is read in blocks, whose coordinates
  do not separate;

- a Laplace or elastic-net parent not centered at zero, the closed form
  being written for one that is;

- SCAD at `step >= a - 1` and MCP at `step >= gamma`, where the
  subproblem is not convex and the operator is set-valued;

- a `step` that is not a single positive number, or a `v` whose length
  is not `pen@n_coef`.

## Methods

The method on the base class
[`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
signals an error naming the penalty and pointing at
[`has_prox()`](https://statmodels7.github.io/penalties7/reference/has_prox.md).
A penalty supplies the operator when it is quadratic, or separable under
a map that is the identity or diagonal; a separable penalty under a
general map is the generalized-lasso problem and has no closed operator,
so it is rejected rather than approximated.

## See also

[`has_prox()`](https://statmodels7.github.io/penalties7/reference/has_prox.md)
to ask before calling,
[`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
for the same operator as a table a compiled loop can read,
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
for where a gradient method would fail,
[`penalty_gradient()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
for what to use where the penalty is smooth.

## Examples

``` r
v <- c(2, 0.3, -1.4)

# The lasso: a soft threshold, so the middle coordinate becomes exactly 0.
penalty_prox(lasso_penalty(n_coef = 3), v, 1, list(lambda = 1))
#> [1]  1.0  0.0 -0.4

# A ridge: a shrinkage towards zero, and nothing reaches it.
penalty_prox(ridge_penalty(n_coef = 3), v, 1, list(lambda = 1.5))
#> [1]  0.80  0.12 -0.56
v / (1 + 1 * 1.5)
#> [1]  0.80  0.12 -0.56

# The elastic net is the one followed by the other.
penalty_prox(elasticnet_penalty(n_coef = 3), v, 1,
             list(lambda = 1, alpha = 0.6))
#> [1]  1.0000000  0.0000000 -0.5714286

# A heavy-tailed prior has no closed form, so this is a bracketed root.
# It leaves the large coordinate nearly alone, which is the point of it.
penalty_prox(heavy_penalty(n_coef = 3), v, 1, list(sigma = 1, nu = 4))
#> [1]  1.0000000  0.1336635 -0.6579124

# Every one of them minimizes the defining objective. Checking the lasso's
# first coordinate directly:
rho <- lasso_penalty(n_coef = 3)
obj <- function(b1) {
  x <- v; x[1] <- b1
  0.5 * sum((x - v)^2) / 1 + penalty_value(rho, x, list(lambda = 1))
}
stats::optimize(obj, c(-5, 5), tol = 1e-12)$minimum
#> [1] 1

# A branch with no operator says so rather than approximating.
try(penalty_prox(additive_penalty(list(diag(3), diag(c(1, 1, 0)))),
                 v, 1, list(lambda1 = 1, lambda2 = 1)))
#> Error : 'additive [2 components]' has no proximal operator.
#>   A penalty supplies one when it is quadratic, or separable under the
#>   identity map. Use has_prox() to ask before calling.
```
