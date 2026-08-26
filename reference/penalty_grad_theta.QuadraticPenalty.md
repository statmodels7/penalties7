# Hyperparameter Derivatives of a Quadratic Penalty

The three blocks in \\\lambda\\, all closed form.
[`penalty_grad_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
returns \\\tfrac{1}{2}(D\beta)'P(D\beta) - r/(2\lambda)\\,
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
returns \\r/(2\lambda^2)\\, and
[`penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
returns \\D'PD\beta\\, which is the gradient with \\\lambda\\ divided
out.

## Arguments

- pen:

  A
  [`QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/QuadraticPenalty.md)
  object.

- beta:

  A numeric vector of length `pen@n_coef`.

- theta:

  A list holding `lambda`.

- scale:

  Read by the generic, which applies the chain rule onto the
  unconstrained scale after dispatch. These methods always return the
  parameter scale.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`penalty_grad_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
a list of one number named `lambda`;
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
a list of one number named `lambda_lambda`;
[`penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
a list of one numeric vector of length `pen@n_coef`, named `lambda`.

## Details

The value is affine in \\\lambda\\ apart from the
\\-\tfrac{r}{2}\log\lambda\\ of the normalizing constant, so the second
derivative comes from that term alone and is positive: the value is
convex in the smoothing parameter, with its minimum at \\\lambda = r /
(D\beta)'P(D\beta)\\. That is the joint-mode estimate of \\\lambda\\ at
fixed coefficients, and it is finite only because the constant was kept.

The mixed block carries no \\\lambda\\ at all, the value being linear in
it in every coefficient direction.

## See also

[`penalty_value.QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.QuadraticPenalty.md)
for the quantity differentiated,
[`penalty_dhessian.QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.QuadraticPenalty.md)
for the third-order blocks,
[`penalty_logpdet.QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.QuadraticPenalty.md)
for the constant's derivatives.

## Examples

``` r
pen <- quadratic_penalty(diag(3))
b <- c(1, -0.5, 0.2)

# The value is convex in lambda, with its minimum where the gradient
# vanishes: at r divided by the quadratic form.
lam_hat <- 3 / sum(b^2)
penalty_grad_theta(pen, b, list(lambda = lam_hat))
#> $lambda
#> [1] 0
#> 
penalty_hess_theta(pen, b, list(lambda = lam_hat))
#> $lambda_lambda
#> [1] 0.27735
#> 

# The mixed block is the gradient with lambda divided out.
all.equal(penalty_cross(pen, b, list(lambda = 4))$lambda,
          penalty_gradient(pen, b, list(lambda = 4)) / 4)
#> [1] TRUE
```
