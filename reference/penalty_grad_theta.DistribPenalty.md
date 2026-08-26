# Hyperparameter Derivatives of a Separable Penalty

The parent's own score and information in its parameters, negated and
summed over the blocks, plus the mixed block.
[`penalty_grad_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
returns \\-\sum_i \ell^{(\theta_k)}\\,
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
returns \\-\sum_i \ell^{(\theta_k\theta_l)}\\, and
[`penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
returns \\-D'\ell^{(y\theta_k)}\\, one coefficient vector per
hyperparameter.

## Arguments

- pen:

  A
  [`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
  object.

- beta:

  A numeric vector of length `pen@n_coef`.

- theta:

  A named list of the parent's free parameters.

- scale:

  Read by the generic, which applies the chain rule onto the
  unconstrained scale after dispatch, using the parent's own links.
  These methods always return the parameter scale.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`penalty_grad_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
a list of one number per hyperparameter, named by `pen@params`.
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
a list of one number per unordered pair, keyed diagonals first.
[`penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
a list of one numeric vector of length `pen@n_coef` per hyperparameter.

## Details

The first two are
[`distributions7::distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.html)
and
[`distributions7::distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.html)
read at the blocks and summed, the map not entering: a hyperparameter of
the prior does not act through \\D\\. The third is
[`distributions7::distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.html),
which exists for this, carried back through the map like the coefficient
gradient.

## The reordering

distributions7 keys its Hessian components lexicographically and this
package keys them diagonals first, so the second derivatives are looked
up under both spellings of each pair, `a_b` and `b_a`, and returned
under this package's. The key is composed from the pair rather than
parsed out of a name, a parameter name being free to contain an
underscore.

## See also

[`penalty_value.DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.DistribPenalty.md)
for the quantity differentiated,
[`penalty_grad_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
for the generic and its two scales,
[`distributions7::distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.html)
for the mixed block's source.

## Examples

``` r
b <- c(1, 0, -1)
pen <- lasso_penalty(n_coef = 3)

# The lasso's value is lambda * L1 - q log(lambda/2), so its derivative in
# lambda is L1 - q/lambda and the second is q/lambda^2.
penalty_grad_theta(pen, b, list(lambda = 0.5))
#> $lambda
#> [1] -4
#> 
sum(abs(b)) - 3 / 0.5
#> [1] -4
penalty_hess_theta(pen, b, list(lambda = 0.5))
#> $lambda_lambda
#> [1] 12
#> 
3 / 0.5^2
#> [1] 12

# The mixed block is the sign vector, the value being lambda times the L1
# norm in every coefficient direction.
penalty_cross(pen, b, list(lambda = 0.5))
#> $lambda
#> [1]  1  0 -1
#> 

# An elastic net has two hyperparameters, so three Hessian keys.
names(penalty_hess_theta(elasticnet_penalty(n_coef = 3), b,
                         list(lambda = 1, alpha = 0.4)))
#> [1] "lambda_lambda" "alpha_alpha"   "lambda_alpha" 
```
