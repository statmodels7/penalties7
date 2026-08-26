# Marginal Derivatives of a Separable Penalty

One page for the branch's answers to the four generics a marginal
criterion asks beyond the second order. With \\\rho = -\sum_j \log
f((D\beta)\_j;\theta)\\ the Hessian is
\\-D'\mathrm{diag}(\ell^{(yy)})D\\, so its \\\theta\\-derivatives are
the parent's own higher components carried through the same map.

## Arguments

- pen:

  A
  [`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
  object.

- beta:

  A numeric vector of length `pen@n_coef`. Read by all three derivative
  methods, the parent's components depending on the argument.

- theta:

  A named list of the parent's free parameters.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
a list of one matrix of side `pen@n_coef` per hyperparameter.
[`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md)
a list of one such matrix per unordered pair, keyed diagonals first.
[`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md)
a list of one numeric vector of length `pen@n_coef` per unordered pair.
[`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md)
a single logical.

## Where each comes from

|  |  |
|----|----|
| generic | parent's component |
| [`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md) | [`distributions7::distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.html) |
| [`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md) | [`distributions7::distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.html) |
| [`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md) | [`distributions7::distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.html) |

Nothing is differentiated here. A parent with closed forms for those,
the gaussian among them and so every Gaussian random effect, makes this
branch exact; a parent without them inherits that package's documented
fallback, one central difference of its analytic first-order component.

A multivariate parent is read blockwise: each component is assembled
into a block-diagonal matrix by
[`dp_blockdiag()`](https://statmodels7.github.io/penalties7/reference/dp_blockdiag.md)
and carried through the map by
[`map_quad_full()`](https://statmodels7.github.io/penalties7/reference/map_quad.md),
where a univariate one goes through
[`map_quad()`](https://statmodels7.github.io/penalties7/reference/map_quad.md).

## What rejects

A penalty whose parent declares a kink, which is the lasso and the
elastic net. The third derivative does not exist at the kink, and the
mode a marginal criterion expands around is exactly where a kinked
penalty puts coefficients, so all three derivative generics reject with
the penalty named.
[`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md)
does not: it answers `TRUE` for the lasso, whose log-density is linear
in the response away from the kink.

## How [`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md) decides

By probing the parent's
[`distributions7::distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.html)
at four points and asking whether it is constant to `1e-12`. The
log-density is quadratic in the response exactly when its second
derivative there does not depend on it, and that second derivative is
analytic for almost every family, where the third is often a difference
whose noise no threshold separates from a true zero.

## See also

[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
and its siblings for the generics,
[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md)
for the branch,
[`reject_kinked()`](https://statmodels7.github.io/penalties7/reference/reject_kinked.md)
for the refusal,
[`carry_pairs()`](https://statmodels7.github.io/penalties7/reference/carry_pairs.md)
for the re-keying.

## Examples

``` r
b <- c(1, -0.5, 0.3)

# A Gaussian prior at zero: the Hessian is 1/sigma^2 on the diagonal, so
# its derivative in sigma is -2/sigma^3.
g <- distrib_penalty(
  distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0),
  n_coef = 3)
diag(penalty_dhessian(g, b, list(sigma = 1.2))$sigma)
#> [1] -1.157407 -1.157407 -1.157407
rep(-2 / 1.2^3, 3)
#> [1] -1.157407 -1.157407 -1.157407

# A Student t prior is not quadratic in the coefficients, so its
# derivative moves with them.
h <- heavy_penalty(n_coef = 3)
beta_quadratic(h, list(sigma = 1, nu = 4))
#> [1] FALSE
diag(penalty_dhessian(h, b, list(sigma = 1, nu = 4))$nu)
#> [1]  0.080000000 -0.004070832 -0.038850928

# A kinked parent has no third derivative and says so.
try(penalty_dhessian(lasso_penalty(n_coef = 3), b, list(lambda = 1)))
#> Error : 'separable [fixed laplace2 [mu=0]]' has a kink, so penalty_dhessian() does not exist there and
#>   its hyperparameters cannot be estimated by a marginal criterion.
```
