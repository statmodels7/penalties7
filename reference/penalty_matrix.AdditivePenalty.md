# Marginal Quantities of an Additive Penalty

The four quantities that a marginal criterion reads.
[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
returns the weighted sum \\S(\lambda) = \sum_k \lambda_k P_k\\,
[`penalty_rank()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
the rank and
[`penalty_null_basis()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
an orthonormal basis of the null space shared by the components, both
fixed at construction, and
[`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
the log pseudo-determinant with its first two derivatives.

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
an orthonormal basis of the null space shared by the components, a
matrix with `pen@n_coef` rows and `pen@n_coef - penalty_rank(pen)`
columns.
[`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
a list of `value` (a single number), `grad` (a list keyed by
`pen@params`) and `hess` (a list keyed by the pairs that
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
uses). Every entry of it is `NaN` at a setting that
[`additive_sum()`](https://statmodels7.github.io/penalties7/reference/additive_sum.md)
reports as not resolvable.

## Details

With \\S^{+}\\ the pseudo-inverse over the stored rank,

\$\$\frac{\partial}{\partial\lambda_k}\log\mathrm{pdet}\\S =
\operatorname{tr}(S^{+}P_k), \qquad
\frac{\partial^{2}}{\partial\lambda_k\partial\lambda_l}
\log\mathrm{pdet}\\S = -\operatorname{tr}(S^{+}P_kS^{+}P_l).\$\$

Both formulas are exact. The value, the gradient and the Hessian of the
log pseudo-determinant come from the similarity transformation of
[`additive_sum()`](https://statmodels7.github.io/penalties7/reference/additive_sum.md),
which takes the traces in the transformed coordinates and not from a
materialized \\S^{+}\\. All three are `NaN` at a setting that
[`additive_sum()`](https://statmodels7.github.io/penalties7/reference/additive_sum.md)
reports as not resolvable.
[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
assembles the sum directly, and
[`penalty_rank()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
and
[`penalty_null_basis()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
return the stored values.

[`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
returns the same shape as
[`penalty_logpdet.QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.QuadraticPenalty.md)
and
[`penalty_logpdet.StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.StructuredPenalty.md):
`grad` is a list keyed by `pen@params` and `hess` a list keyed by the
pairs that
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
uses, so code written for those branches reads `grad$lambda1` here.

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
