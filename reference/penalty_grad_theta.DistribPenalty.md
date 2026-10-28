# Hyperparameter Derivatives of a Separable Penalty

These methods return the parent's own score and information in its
parameters, negated and summed over the blocks, plus the mixed block.
[`penalty_grad_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
returns \\-\sum_i \ell^{(\theta_k)}\\,
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
returns \\-\sum_i \ell^{(\theta_k\theta_l)}\\, and
[`penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
returns \\-D'\ell^{(y\theta_k)}\\, one coefficient vector per
hyperparameter.

## Arguments

- penalty:

  A
  [`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
  object.

- beta:

  A numeric vector of length `penalty@n_coef`.

- theta:

  A named list of the parent's free parameters.

- scale:

  Read by the generic, which applies the chain rule to the unconstrained
  scale after dispatch and uses the links copied from the parent. These
  methods always return derivatives on the parameter scale.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`penalty_grad_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
returns a list of one number per hyperparameter, named by
`penalty@params`.
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
returns a list of one number per unordered pair, keyed diagonals first.
[`penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
returns a list of one numeric vector of length `penalty@n_coef` per
hyperparameter.

## Details

The first two are
[`distributions7::distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.html)
and
[`distributions7::distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.html)
evaluated at the blocks and summed. The map does not enter them, because
a hyperparameter of the prior does not act through \\D\\. The third is
[`distributions7::distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.html),
carried back through the map like the coefficient gradient.

## Keys of the second derivatives

The second derivatives are returned under the keys of
[`ptheta_pairs()`](https://statmodels7.github.io/penalties7/reference/ptheta_pairs.md):
the diagonal pairs first, then the off-diagonal pairs, in the order of
`penalty@params`. The shipped continuous families of distributions7 use
the same order, but a family may spell a pair the other way, so each
pair is looked up under both spellings, `a_b` and `b_a`. The key is
composed from the pair and not parsed out of a name, because a parameter
name may contain an underscore.

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

# The mixed block is the sign vector, because the value is lambda times the
# L1 norm in every coefficient direction.
penalty_cross(pen, b, list(lambda = 0.5))
#> $lambda
#> [1]  1  0 -1
#> 

# An elastic net has two hyperparameters, so three Hessian keys.
names(penalty_hess_theta(elasticnet_penalty(n_coef = 3), b,
                         list(lambda = 1, alpha = 0.4)))
#> [1] "lambda_lambda" "alpha_alpha"   "lambda_alpha" 
```
