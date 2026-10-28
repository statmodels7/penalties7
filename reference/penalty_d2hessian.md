# The Second Derivative of the Coefficient Hessian in the Hyperparameters

Returns \\\partial^2 S/\partial\theta_m\partial\theta_l =
\partial^4\rho/\partial\beta^2\partial\theta_m\partial\theta_l\\, one
matrix per unordered pair. The second derivative of a marginal criterion
in the hyperparameters contains it, and so does an exact outer Hessian.

## Usage

``` r
penalty_d2hessian(penalty, beta, theta, ...)
```

## Arguments

- penalty:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- beta:

  A numeric vector of length `penalty@n_coef`; the generic signals an
  error for any other length. Only the separable branch reads its
  values.

- theta:

  A named list of hyperparameter values, or a named numeric vector
  carrying the same, checked by the generic as in
  [`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md).

- ...:

  Passed to methods. The shipped methods do not read it.

## Value

A named list of one matrix per unordered pair of hyperparameters, each
of side `penalty@n_coef`, keyed diagonals first. A penalty with \\p\\
hyperparameters gives \\p(p+1)/2\\ entries.

## Details

It is exactly zero for every penalty whose Hessian is linear in its
hyperparameters, which is the quadratic branch, where \\S = \lambda
D'PD\\, and the additive branch, where \\S = \sum_k \lambda_k D'P_kD\\.
The structured branch reads the second derivative `param_d2` of the
matrix parameter. The separable branch reads the second
\\\theta\\-derivative of the response curvature of the parent through
[`distributions7::distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.html)
and carries it through the map as in
[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md).
The derivatives are taken in \\\theta\\ on its own scale.

The keys are those of
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md):
the diagonal pairs come first in `penalty@params` order, then the upper
off-diagonal pairs follow, and each key joins the two names with an
underscore. Each unordered pair has one key, written in `penalty@params`
order (`sigma_nu` and not `nu_sigma` when `sigma` comes first), and the
reversed name is not an element of the list.

A penalty with a kink signals an error, and so does the base class.

## See also

[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
for the first order,
[`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md)
for the mixed one,
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
for the keying.

## Examples

``` r
b <- c(1, -0.5, 0.3)

# The second derivative is zero for a quadratic penalty, whose Hessian is
# linear in lambda.
pen <- quadratic_penalty(diag(3))
penalty_d2hessian(pen, b, list(lambda = 2))
#> $lambda_lambda
#>      [,1] [,2] [,3]
#> [1,]    0    0    0
#> [2,]    0    0    0
#> [3,]    0    0    0
#> 

# It is zero for an additive penalty as well, and the keys list the
# diagonals first over three pairs.
add <- additive_penalty(list(diag(3), diag(c(1, 1, 0))))
names(penalty_d2hessian(add, b, list(lambda1 = 2, lambda2 = 0.5)))
#> [1] "lambda1_lambda1" "lambda2_lambda2" "lambda1_lambda2"

# It is not zero for a structured penalty, whose matrix is not linear in
# its free values.
s <- structured_penalty(parameters7::log_cholesky(2))
th <- list(log_L1 = 0.1, log_L2 = -0.1, L2.1 = 0.3)
penalty_d2hessian(s, c(1, -0.5), th)$log_L1_log_L1
#>           [,1]      [,2]
#> [1,] 4.8856110 0.3315513
#> [2,] 0.3315513 0.0000000
```
