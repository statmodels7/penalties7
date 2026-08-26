# Coefficient Derivatives of an Additive Penalty

[`penalty_gradient()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
returns \\S(\lambda)\beta\\ and
[`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
returns \\S(\lambda)\\. Both assemble the weighted sum directly and take
no eigendecomposition, the normalizing constant not depending on the
coefficients.

## Arguments

- pen:

  An
  [`AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/AdditivePenalty.md)
  object.

- beta:

  A numeric vector of length `pen@n_coef`.

- theta:

  A named list of `lambda1`, `lambda2`, ...

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`penalty_gradient()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
a numeric vector of length `pen@n_coef`;
[`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
a symmetric base matrix of that side.

## Details

The value is a quadratic form, so the Hessian does not move with the
coefficients and the gradient is the Hessian applied to them. A
coefficient vector lying in the intersected null space of the components
has zero gradient at every setting of the parameters.

## See also

[`penalty_value.AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.AdditivePenalty.md)
for the quantity differentiated,
[`penalty_grad_theta.AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.AdditivePenalty.md)
for the hyperparameter blocks.

## Examples

``` r
P1 <- crossprod(diff(diag(5)))
P2 <- crossprod(diff(diag(5), differences = 2))
pen <- additive_penalty(list(P1, P2))
th <- list(lambda1 = 2, lambda2 = 0.5)
b <- c(0.4, -1.1, 0.7, 0.2, -0.3)

all.equal(penalty_gradient(pen, b, th),
          drop(penalty_hessian(pen, b, th) %*% b))
#> [1] TRUE

# The constants are in both null spaces, so they cost nothing at any
# setting of the two parameters.
penalty_gradient(pen, rep(1, 5), th)
#> [1] 0 0 0 0 0
penalty_gradient(pen, rep(1, 5), list(lambda1 = 1e6, lambda2 = 1e-6))
#> [1] 0 0 0 0 0
```
