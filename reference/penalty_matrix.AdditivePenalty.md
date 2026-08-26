# Marginal Quantities of an Additive Penalty

Three of the four pieces a marginal criterion reads.
[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
returns the weighted sum \\S(\lambda) = \sum_k \lambda_k P_k\\,
[`penalty_rank()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
the rank fixed at construction, and
[`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
the log pseudo-determinant with its first two derivatives. There is no
[`penalty_null_basis()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
method for this branch, and the base class's rejects.

## Arguments

- pen:

  An
  [`AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/AdditivePenalty.md)
  object.

- theta:

  A named list of `lambda1`, `lambda2`, ...
  [`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
  and
  [`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
  only.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
a symmetric base matrix of side `pen@n_coef`.
[`penalty_rank()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
a single integer, and
[`penalty_null_basis()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
an orthonormal basis of the components' shared null space, `pen@n_coef`
by `pen@n_coef - penalty_rank(pen)`.
[`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
a list of `value` (a single number), `grad` (a list keyed by
`pen@params`) and `hess` (a list keyed by the pairs
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
uses), which is the shape the quadratic and structured branches answer
in.

## Details

With \\S^{+}\\ the pseudo-inverse over the stored rank,

\$\$\frac{\partial}{\partial\lambda_k}\log\mathrm{pdet}\\S =
\operatorname{tr}(S^{+}P_k), \qquad
\frac{\partial^{2}}{\partial\lambda_k\partial\lambda_l}
\log\mathrm{pdet}\\S = -\operatorname{tr}(S^{+}P_kS^{+}P_l),\$\$

both exact and both agreeing with the traces computed apart to 0.

**[`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
answers in a different shape here.** `grad` is an unnamed numeric vector
in `pen@params` order and `hess` is a square matrix, where
[`penalty_logpdet.QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.QuadraticPenalty.md)
and
[`penalty_logpdet.StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.StructuredPenalty.md)
return named lists keyed by hyperparameter and by pair. A consumer
written against those will read `NULL` from `grad$lambda1` here. Nothing
in the toolkit reads it today, because
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
answers `FALSE` for this branch and every consumer routes on that first.

## See also

[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
and its siblings for the generics,
[`additive_sum()`](https://statmodels7.github.io/penalties7/reference/additive_sum.md)
for the decomposition these read,
[`penalty_matrix.QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.QuadraticPenalty.md)
for the branch whose determinant is linear in one parameter.

## Examples

``` r
P1 <- crossprod(diff(diag(5)))
P2 <- crossprod(diff(diag(5), differences = 2))
pen <- additive_penalty(list(P1, P2))
th <- list(lambda1 = 2, lambda2 = 0.5)

penalty_rank(pen)
#> [1] 4
lp <- penalty_logpdet(pen, th)
str(lp)
#> List of 3
#>  $ value: num 5.54
#>  $ grad :List of 2
#>   ..$ lambda1: num 1.54
#>   ..$ lambda2: num 1.84
#>  $ hess :List of 3
#>   ..$ lambda1_lambda1: num -0.624
#>   ..$ lambda2_lambda2: num -1.36
#>   ..$ lambda1_lambda2: num -0.582

# The gradient is the trace of the pseudo-inverse against each component.
S <- penalty_matrix(pen, th)
e <- eigen(S, symmetric = TRUE)
k <- order(e$values, decreasing = TRUE)[seq_len(penalty_rank(pen))]
Sp <- e$vectors[, k, drop = FALSE] %*% (t(e$vectors[, k, drop = FALSE]) /
                                          e$values[k])
max(abs(unlist(lp$grad) - c(sum(Sp * P1), sum(Sp * P2))))
#> [1] 0

# The null space of the sum is the intersection of the components', so it
# does not move with the hyperparameters.
nb <- penalty_null_basis(pen)
c(width = ncol(nb), annihilated = max(abs(S %*% nb)) < 1e-12)
#>       width annihilated 
#>           1           1 

# And the value is the sum of the logarithms of the non-zero eigenvalues.
lp$value - sum(log(e$values[k]))
#> [1] 0
```
