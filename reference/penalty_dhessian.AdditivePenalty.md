# Marginal Derivatives of an Additive Penalty

The methods of the additive branch cover the four generics that a
marginal criterion requires beyond the second order. With \\S_k =
D'P_kD\\, the Hessian \\S = \sum_k \lambda_k S_k\\ is linear in the
smoothing parameters and free of the coefficients, so
[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
returns the components \\S_k\\ themselves,
[`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md)
and
[`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md)
return zeros, and
[`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md)
returns `TRUE`.

## Arguments

- pen:

  An
  [`AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/AdditivePenalty.md)
  object.

- beta:

  A numeric vector of length `pen@n_coef`. The methods do not use it.

- theta:

  A named list of `lambda1`, `lambda2`, ... The methods do not use it.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
returns a list of one base matrix per component, named by `pen@params`.
[`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md)
returns a list of zero matrices, one per unordered pair.
[`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md)
returns a list of zero vectors, one per unordered pair.
[`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md)
returns the single logical `TRUE`.

## Details

Every entry is exact, and the methods do not compute any entry at call
time. The components were fixed at construction and are stored already
carried through the map, so the first derivative is a lookup of
`pen@mats` and the higher ones are zeros of the right shape. Under the
identity map \\S_k = P_k\\.

[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
returns `TRUE` for this branch, as it does for the quadratic and
structured branches.

## See also

[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
and its siblings for the generics,
[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
for the branch,
[`penalty_dhessian_beta.AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian_beta.penalty.md)
and
[`penalty_d2hessian_beta.AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.penalty.md)
for the derivatives of the Hessian in the coefficients.

## Examples

``` r
add <- additive_penalty(list(crossprod(diff(diag(4))), diag(4)))
b <- c(1, -0.5, 0.3, 0.2)
th <- list(lambda1 = 2, lambda2 = 0.5)

# The derivative in each parameter is that parameter's own component.
d <- penalty_dhessian(add, b, th)
names(d)
#> [1] "lambda1" "lambda2"
max(abs(d$lambda2 - diag(4)))
#> [1] 0

# Everything above first order is zero, over all three pairs.
names(penalty_d2hessian(add, b, th))
#> [1] "lambda1_lambda1" "lambda2_lambda2" "lambda1_lambda2"
max(abs(unlist(penalty_d2hessian(add, b, th))))
#> [1] 0
```
