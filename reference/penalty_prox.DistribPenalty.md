# Proximal Operator of a Separable Penalty

Closed form for the three families that have one, and a coordinatewise
root for every other parent with a differentiable log-density. A
diagonal map is handled by rescaling; anything else is rejected.

## Arguments

- pen:

  A
  [`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
  object.

- v:

  A numeric vector of length `pen@n_coef`.

- step:

  The step length \\t\\, a single positive number.

- theta:

  A named list of hyperparameter values.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of the same length as `v`.

## The closed forms

|  |  |
|----|----|
| parent | operator |
| `gaussian1` | \\(v - t g_0)/(1 + t/\sigma^2)\\, with \\g_0\\ the gradient at the origin |
| `laplace`, `laplace2` | \\\mathrm{sign}(v)(\lvert v\rvert - t\lambda)\_{+}\\ |
| `enet` | the soft threshold at \\t\lambda\alpha\\, then division by \\1 + t\lambda(1-\alpha)\\ |

The elastic net's two steps come out of one stationary condition:
\\(\beta - v)/t + a\\\mathrm{sign}(\beta) + c\beta = 0\\ with \\a =
\lambda\alpha\\ and \\c = \lambda(1-\alpha)\\ separates into the Laplace
part's threshold followed by the Gaussian part's shrinkage.

## Where the parent must sit

The Laplace and elastic-net forms are written for a parent centered at
zero, and a family name does not say where a location is. The gradient
at the origin is evaluated and required to vanish, and an off-center
parent is rejected. The Gaussian needs no such restriction: its location
enters the stationary condition linearly and is carried through as
\\g_0\\.

## The root

Any other parent is solved coordinatewise from \\(\beta - v)/t =
\ell^{(y)}(\beta)\\, the right-hand side being
[`distributions7::distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.html).
For a log-concave density the difference of the two sides is strictly
increasing, so the root is unique; the bracket starts one unit either
side of \\\min(v, 0)\\ and \\\max(v, 0)\\ and is doubled up to sixty
times, then [`uniroot()`](https://rdrr.io/r/stats/uniroot.html) refines
to `eps^0.75`. A parent declaring a kink, and outside the closed-form
families, is rejected instead: the root is a smooth instrument and the
kink is where the answer is.

## Errors

A parent read in blocks, a map that is not diagonal, an off-center
Laplace or elastic-net parent, a parent with an unrecognized kink, and a
bracket that could not be found, each with a message naming the penalty
and the reason.

## See also

[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
for the generic,
[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md)
for the branch,
[`distributions7::distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.html)
for the root's right-hand side.

## Examples

``` r
v <- c(2, 0.3, -1.4)

# The soft threshold sets everything within t * lambda of zero to zero.
penalty_prox(lasso_penalty(n_coef = 3), v, 1, list(lambda = 1))
#> [1]  1.0  0.0 -0.4

# The elastic net thresholds first and then shrinks what survives.
penalty_prox(elasticnet_penalty(n_coef = 3), v, 1,
             list(lambda = 1, alpha = 0.6))
#> [1]  1.0000000  0.0000000 -0.5714286
sign(v) * pmax(abs(v) - 1 * 1 * 0.6, 0) / (1 + 1 * 1 * (1 - 0.6))
#> [1]  1.0000000  0.0000000 -0.5714286

# A Student t prior has no closed form. The root leaves the large
# coordinate nearly where it was, which a Gaussian prior would not.
penalty_prox(heavy_penalty(n_coef = 3), v, 1, list(sigma = 1, nu = 4))
#> [1]  1.0000000  0.1336635 -0.6579124
penalty_prox(ridge_penalty(n_coef = 3), v, 1, list(lambda = 1))
#> [1]  1.00  0.15 -0.70
```
