# Coefficient Derivatives of a Quadratic Penalty

[`penalty_gradient()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
returns \\\lambda D'PD\beta\\ and
[`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
returns \\\lambda D'PD\\. Both read the cached \\D'PD\\, so neither
touches the map at call time, and the Hessian does not depend on
\\\beta\\ at all, the value being a quadratic form.

## Arguments

- pen:

  A
  [`QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/QuadraticPenalty.md)
  object.

- beta:

  A numeric vector of length `pen@n_coef`.

- theta:

  A list holding `lambda`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`penalty_gradient()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
a numeric vector of length `pen@n_coef`;
[`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
a symmetric matrix of that side, a base matrix for an ordinary penalty
and a `dgCMatrix` when the penalty was built with `blocks > 1`.

## Details

Because the Hessian is constant, the gradient is exactly the Hessian
applied to the coefficients, and a Newton step on this penalty alone
reaches the minimum in one iteration from anywhere.

## See also

[`penalty_value.QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.QuadraticPenalty.md)
for the quantity differentiated,
[`penalty_grad_theta.QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.QuadraticPenalty.md)
for the hyperparameter derivatives.

## Examples

``` r
pen <- quadratic_penalty(crossprod(diff(diag(4))))
b <- c(1, -0.5, 0.2, 2)

# The gradient is the constant Hessian applied to the coefficients.
all.equal(penalty_gradient(pen, b, list(lambda = 3)),
          drop(penalty_hessian(pen, b, list(lambda = 3)) %*% b))
#> [1] TRUE

# The Hessian does not move with beta.
all.equal(penalty_hessian(pen, b, list(lambda = 3)),
          penalty_hessian(pen, 0 * b, list(lambda = 3)))
#> [1] TRUE
```
