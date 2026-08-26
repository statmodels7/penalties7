# Marginal Derivatives of an Additive Penalty

One page for the branch's answers to the four generics a marginal
criterion asks beyond the second order. \\S = \sum_k \lambda_k P_k\\ is
linear in the smoothing parameters and free of the coefficients, so
[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
returns the components themselves,
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

  A numeric vector of length `pen@n_coef`. Unused.

- theta:

  A named list of `lambda1`, `lambda2`, ... Unused.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
a list of one matrix per component, named by `pen@params`.
[`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md)
a list of zero matrices, one per unordered pair.
[`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md)
a list of zero vectors, one per unordered pair.
[`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md)
the single logical `TRUE`.

## Details

Every entry is exact and none is computed at call time: the components
were fixed at construction, already carried through the map, so the
first derivative is a lookup and the higher ones are zeros of the right
shape.

This is the one branch whose marginal derivatives are reachable while
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
answers `FALSE` for it. A consumer routing on that predicate rather than
on these methods will not find them; see
[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md).

## See also

[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
and its siblings for the generics,
[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
for the branch and for what
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
costs it.

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
