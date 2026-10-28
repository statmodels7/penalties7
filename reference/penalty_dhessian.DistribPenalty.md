# Marginal Derivatives of a Separable Penalty

The methods of the separable branch cover the four generics that a
marginal criterion requires beyond the second order. With \\\rho =
-\sum_j \log f((D\beta)\_j;\theta)\\ the Hessian is
\\-D'\mathrm{diag}(\ell^{(yy)})D\\, so its \\\theta\\-derivatives are
the higher components of the parent carried through the same map.

## Arguments

- penalty:

  A
  [`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
  object.

- beta:

  A numeric vector of length `penalty@n_coef`. All three derivative
  methods read it, because the components of the parent depend on it.

- theta:

  A named list of the free parameters of the parent.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
returns a list of one base matrix of side `penalty@n_coef` per
hyperparameter.
[`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md)
returns a list of one such matrix per unordered pair, keyed diagonals
first.
[`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md)
returns a list of one numeric vector of length `penalty@n_coef` per
unordered pair.
[`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md)
returns a single logical.

## The component of the parent behind each method

|  |  |  |
|----|----|----|
| generic | component of the parent | result |
| [`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md) | [`distributions7::distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.html) | \\-D'\mathrm{diag}(\partial^3\ell/\partial y^2\partial\theta_m)D\\ |
| [`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md) | [`distributions7::distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.html) | \\-D'\mathrm{diag}(\partial^4\ell/\partial y^2\partial\theta_m\partial\theta_l)D\\ |
| [`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md) | [`distributions7::distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.html) | \\-D'(\partial^3\ell/\partial y\\\partial\theta_m\partial\theta_l)\\ |

Each component is evaluated at \\D\beta\\ and at `theta` on the
parameter scale of the parent, so the derivatives are in \\\theta\\
itself, without a link.

This branch does not differentiate. A parent with closed forms for these
components, such as the Gaussian and therefore every Gaussian random
effect, makes this branch exact. A parent without them takes the
fallback of distributions7, which applies one stencil to the
highest-order quantity that the family implements itself, usually its
analytic response Hessian.

A multivariate parent is read blockwise: each component is assembled
into a block-diagonal matrix by
[`dp_blockdiag()`](https://statmodels7.github.io/penalties7/reference/dp_blockdiag.md)
and carried through the map by
[`map_quad_full()`](https://statmodels7.github.io/penalties7/reference/map_quad.md),
and a univariate one goes through
[`map_quad()`](https://statmodels7.github.io/penalties7/reference/map_quad.md).

## Kinked parents

A penalty whose parent declares a kink, which is the lasso and the
elastic net, signals an error in
[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md),
[`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md)
and
[`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md),
with the penalty named (see
[`reject_kinked()`](https://statmodels7.github.io/penalties7/reference/reject_kinked.md)).
The third derivative does not exist at the kink, and a kinked penalty
places coefficients at the mode at which a marginal criterion is
expanded.
[`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md)
does not signal an error and returns `TRUE` for the lasso, whose
log-density is linear in the response away from the kink.

## The test of [`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md)

The log-density is quadratic in the response exactly when its second
derivative in the response does not depend on the response. The generic
has already checked `theta` (an invalid value is an error), and the
method evaluates
[`distributions7::distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.html)
of the parent at four points of `penalty@block` coordinates each, and
tests with [`all.equal()`](https://rdrr.io/r/base/all.equal.html) at a
relative tolerance of `1e-12` whether the values agree across the
points. A multivariate parent returns one \\p \times p\\ matrix per
point, or a single matrix when that matrix does not depend on the
response, and the method compares the matrices of different points with
one another, not the entries of one matrix. A single value or a single
matrix, which carries no dimension for the points, gives `TRUE`. If the
parent signals an error, or returns a number of values that is not a
multiple of `penalty@block^2`, the method returns `FALSE`.

## See also

[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
and its siblings for the generics,
[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md)
for the branch,
[`reject_kinked()`](https://statmodels7.github.io/penalties7/reference/reject_kinked.md)
for the error on a kinked parent,
[`carry_pairs()`](https://statmodels7.github.io/penalties7/reference/carry_pairs.md)
for the re-keying,
[`penalty_dhessian_beta.DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian_beta.penalty.md)
and
[`penalty_d2hessian_beta.DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.penalty.md)
for the derivatives of the Hessian in the coefficients.

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
# derivative changes with them.
h <- heavy_penalty(n_coef = 3)
beta_quadratic(h, list(sigma = 1, nu = 4))
#> [1] FALSE
diag(penalty_dhessian(h, b, list(sigma = 1, nu = 4))$nu)
#> [1]  0.080000000 -0.004070832 -0.038850928

# A kinked parent has no third derivative, and the call is an error.
try(penalty_dhessian(lasso_penalty(n_coef = 3), b, list(lambda = 1)))
#> Error : 'separable [fixed laplace2 [mu=0]]' has a kink, so penalty_dhessian() does not exist there and
#>   its hyperparameters cannot be estimated by a marginal criterion.
```
