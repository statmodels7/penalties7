# The Second Derivative of the Coefficient Hessian in the Hyperparameters

Returns \\\partial^2 S/\partial\theta_m\partial\theta_l =
\partial^4\rho/\partial\beta^2\partial\theta_m\partial\theta_l\\, one
matrix per unordered pair. A marginal criterion differentiated twice in
the hyperparameters needs it, and so does an exact outer Hessian.

## Usage

``` r
penalty_d2hessian(pen, beta, theta, ...)
```

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- beta:

  A numeric vector of length `pen@n_coef`. Read only by the separable
  branch.

- theta:

  A named list of hyperparameter values, or a named numeric vector
  carrying the same.

- ...:

  Passed to methods. No shipped method reads it.

## Value

A named list of one matrix per unordered pair of hyperparameters, each
of side `pen@n_coef`, keyed diagonals first. A penalty with \\p\\
hyperparameters gives \\p(p+1)/2\\ entries.

## Details

It is **exactly zero** for every penalty whose Hessian is linear in its
hyperparameters, which is the quadratic branch, where \\S = \lambda
D'PD\\, and the additive branch, where \\S = \sum_k \lambda_k P_k\\. The
structured branch reads the matrix parameter's `param_d2`, and the
separable branch the second \\\theta\\-derivative of the parent's
response curvature, through
[`distributions7::distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.html).

The keys are
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)'s:
diagonals first in `pen@params` order, then the upper off-diagonal
pairs, each joined by an underscore. A consumer looking a pair up
therefore need not know which order it was written in.

A penalty with a kink rejects, and so does the base class.

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

# Zero for a quadratic penalty, whose Hessian is linear in lambda.
pen <- quadratic_penalty(diag(3))
penalty_d2hessian(pen, b, list(lambda = 2))
#> $lambda_lambda
#>      [,1] [,2] [,3]
#> [1,]    0    0    0
#> [2,]    0    0    0
#> [3,]    0    0    0
#> 

# And for an additive one, keyed diagonals first over three pairs.
add <- additive_penalty(list(diag(3), diag(c(1, 1, 0))))
names(penalty_d2hessian(add, b, list(lambda1 = 2, lambda2 = 0.5)))
#> [1] "lambda1_lambda1" "lambda2_lambda2" "lambda1_lambda2"

# Not zero for a structured penalty, whose matrix is not linear in any
# one free value.
s <- structured_penalty(parameters7::log_cholesky(2))
th <- list(log_L1 = 0.1, log_L2 = -0.1, L2.1 = 0.3)
penalty_d2hessian(s, c(1, -0.5), th)$log_L1_log_L1
#>           [,1]      [,2]
#> [1,] 4.8856110 0.3315513
#> [2,] 0.3315513 0.0000000
```
