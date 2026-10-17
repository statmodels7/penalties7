# Proximal Operator of a Separable Penalty

Computes closed forms for the Gaussian, the Laplace (in both
parametrizations) and the elastic net, and for every other parent with a
differentiable log-density it computes the stationary point of smallest
objective, coordinate by coordinate. A diagonal map is handled by
rescaling, and any other map is rejected.

## Arguments

- pen:

  A
  [`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
  object.

- v:

  A numeric vector of length `pen@n_coef`.

- step:

  The step length \\t\\, a single positive finite number.

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

For the `laplace` parametrization \\\lambda = 1/\sigma\\. The two steps
of the elastic net follow from one stationary condition: \\(\beta -
v)/t + a\\\mathrm{sign}(\beta) + c\beta = 0\\, with \\a =
\lambda\alpha\\ and \\c = \lambda(1-\alpha)\\, separates into the
threshold of the Laplace part followed by the shrinkage of the Gaussian
part.

## The location of the parent

The Laplace and elastic-net forms are written for a parent centered at
zero, and the family name does not carry the location. The gradient of
the penalty at the origin is therefore evaluated, and the call signals
an error when its largest absolute entry exceeds \\10^{-8}\max(1,
\lambda)\\. The Gaussian form has no such restriction, because the
location enters the stationary condition linearly and is carried through
as \\g_0\\.

## Parents without a closed form

Any other parent is solved coordinate by coordinate from \\(\beta - v)/t
= \ell^{(y)}(\beta)\\, whose right-hand side is
[`distributions7::distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.html).
The condition can have more than one root when the density is not
log-concave, so every change of sign on the segment between \\v\\ and
the mode of the parent is refined and the root of smallest objective is
returned, as
[`.prox_root()`](https://statmodels7.github.io/penalties7/reference/dot-prox_root.md)
describes. That root is the global minimizer when the parent is
unimodal, and a stationary point that need not be the global minimizer
when the parent has more than one mode. A parent that declares a kink
and is not one of the closed-form families is rejected, because the
stationary condition does not hold at a kink and the minimizer of such a
penalty can lie there.

## Errors

A parent read in blocks, a map that is not a Matrix diagonal object with
nonzero entries, a Laplace or elastic-net parent that is not centered at
zero, a parent that declares a kink and has no closed form, and a
stationary condition with no change of sign on the segment or in the
widened bracket each signal an error with a message giving the reason;
all but the off-center parent name the penalty.

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

# A Student t prior has no closed form, and the operator is found from the
# stationary condition. With sigma = 1, nu = 4 and step 1 it maps 2 to 1,
# as the ridge with lambda = 1 does; the two differ for coordinates that
# are large relative to sigma, such as 10.
v2 <- c(2, 0.3, 10)
penalty_prox(heavy_penalty(n_coef = 3), v2, 1, list(sigma = 1, nu = 4))
#> [1] 1.0000000 0.1336635 9.4958181
penalty_prox(ridge_penalty(n_coef = 3), v2, 1, list(lambda = 1))
#> [1] 1.00 0.15 5.00
```
