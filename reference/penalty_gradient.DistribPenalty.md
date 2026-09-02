# Coefficient Derivatives of a Separable Penalty

The parent's response derivatives, negated and carried back through the
map.
[`penalty_gradient()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
returns \\-D'\ell^{(y)}\\ and
[`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
returns \\-D'\\\mathrm{diag}(\ell^{(yy)})\\D\\, with the middle matrix
block diagonal instead of diagonal when the parent is multivariate.

## Arguments

- pen:

  A
  [`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
  object.

- beta:

  A numeric vector of length `pen@n_coef`.

- theta:

  A named list of the parent's free parameters.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`penalty_gradient()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
a numeric vector of length `pen@n_coef`;
[`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
a symmetric base matrix of that side, diagonal under the identity map
with a univariate parent and block diagonal with a multivariate one.

## Details

Everything comes from
[`distributions7::distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.html)
and
[`distributions7::distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.html),
which are closed form for every continuous family, so no derivative is
taken here.

The `+ 0 * a` in both bodies is a recycling guard. A family whose
derivative does not depend on the argument, a Gaussian's second one for
instance, returns a single number where a vector is wanted, and adding
zero times the argument widens it without changing a value.

At a kink the gradient is one element of the subdifferential, chosen by
the parent. For the lasso at \\\beta_j = 0\\ that is `0`, where the true
subdifferential is \\\[-\lambda, \lambda\]\\;
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
is what produces the exact zero.

## See also

[`penalty_value.DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.DistribPenalty.md)
for the quantity differentiated,
[`penalty_grad_theta.DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.DistribPenalty.md)
for the hyperparameter blocks,
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
for the step that produces exact zeros,
[`dp_blockdiag()`](https://statmodels7.github.io/penalties7/reference/dp_blockdiag.md)
for the multivariate middle matrix.

## Examples

``` r
b <- c(1, 0, -1)

# A Gaussian prior at zero: the gradient is b/sigma^2 and the Hessian is
# constant.
pen <- distrib_penalty(
  distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0),
  n_coef = 3)
penalty_gradient(pen, b, list(sigma = 2))
#> [1]  0.25  0.00 -0.25
b / 4
#> [1]  0.25  0.00 -0.25
diag(penalty_hessian(pen, b, list(sigma = 2)))
#> [1] 0.25 0.25 0.25

# The lasso's gradient is lambda times a sign, and 0 at the kink.
penalty_gradient(lasso_penalty(n_coef = 3), b, list(lambda = 1.5))
#> [1]  1.5  0.0 -1.5

# A multivariate parent gives a block-diagonal Hessian: dependence within
# a block, none between blocks.
mv <- distrib_penalty(
  distributions7::fixed(distributions7::mvgaussian1_distrib(2),
                        mu1 = 0, mu2 = 0), n_coef = 4)
round(penalty_hessian(mv, c(1, 0, -1, 0.5),
                      list(sigma_log_L1 = 0, sigma_log_L2 = 0,
                           sigma_L2.1 = 0.4)), 3)
#>       [,1] [,2]  [,3] [,4]
#> [1,]  1.16 -0.4  0.00  0.0
#> [2,] -0.40  1.0  0.00  0.0
#> [3,]  0.00  0.0  1.16 -0.4
#> [4,]  0.00  0.0 -0.40  1.0
```
