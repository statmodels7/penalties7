# Hyperparameter Derivatives of an Additive Penalty

The three blocks in the smoothing parameters.
[`penalty_grad_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
returns one number per component,
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
one per unordered pair, and
[`penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
one coefficient vector per component.

## Arguments

- pen:

  An
  [`AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/AdditivePenalty.md)
  object.

- beta:

  A numeric vector of length `pen@n_coef`.

- theta:

  A named list of `lambda1`, `lambda2`, ...

- scale:

  Read by the generic, which applies the chain rule onto the
  unconstrained scale after dispatch. These methods always return the
  parameter scale.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`penalty_grad_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
a list of one number per component, named `lambda1`, `lambda2`, ...
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
a list of one number per unordered pair, keyed diagonals first.
[`penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
a list of one numeric vector of length `pen@n_coef` per component.

## Details

Each parameter enters the value linearly through the quadratic form and
non-linearly through the log pseudo-determinant, so with \\S^{+}\\ the
pseudo-inverse of the sum,

\$\$\frac{\partial\rho}{\partial\lambda_k} = \tfrac{1}{2}\beta^\top
P_k\beta - \tfrac{1}{2}\operatorname{tr}(S^{+}P_k), \qquad
\frac{\partial^2\rho}{\partial\lambda_k\partial\lambda_l} =
\tfrac{1}{2}\operatorname{tr}(S^{+}P_kS^{+}P_l), \qquad
\frac{\partial^2\rho}{\partial\beta\\\partial\lambda_k} = P_k\beta.\$\$

The second derivative has no term from the quadratic form, which is
linear in each parameter, and the mixed block has no term from the
determinant.
[`penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
therefore takes no eigendecomposition. The first two blocks read the
traces from
[`additive_sum()`](https://statmodels7.github.io/penalties7/reference/additive_sum.md),
and their entries are `NaN` at a setting that it reports as not
resolvable.

The Hessian is keyed by pairs as
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
keys them: each key joins two parameter names by an underscore
(`lambda1_lambda2`), the diagonal pairs come first and the upper
off-diagonal pairs follow.

## See also

[`penalty_value.AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.AdditivePenalty.md)
for the quantity differentiated,
[`penalty_logpdet.AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.AdditivePenalty.md)
for the determinant's own derivatives,
[`additive_sum()`](https://statmodels7.github.io/penalties7/reference/additive_sum.md)
for the traces that these read.

## Examples

``` r
P1 <- crossprod(diff(diag(5)))
P2 <- crossprod(diff(diag(5), differences = 2))
pen <- additive_penalty(list(P1, P2))
th <- list(lambda1 = 2, lambda2 = 0.5)
b <- c(0.4, -1.1, 0.7, 0.2, -0.3)

penalty_grad_theta(pen, b, th)
#> $lambda1
#> [1] 2.225392
#> 
#> $lambda2
#> [1] 7.168431
#> 
names(penalty_hess_theta(pen, b, th))
#> [1] "lambda1_lambda1" "lambda2_lambda2" "lambda1_lambda2"

# The mixed block is P_k beta, so it does not depend on the parameters.
identical(penalty_cross(pen, b, th),
          penalty_cross(pen, b, list(lambda1 = 99, lambda2 = 1e-3)))
#> [1] TRUE
```
