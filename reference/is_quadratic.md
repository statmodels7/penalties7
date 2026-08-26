# Is a Penalty a Quadratic Form With One Scale?

`TRUE` for
[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
and
[`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md),
the two branches whose value is \\\tfrac{1}{2}(D\beta)^\top S(\theta)
(D\beta)\\ for a matrix that a consumer can ask for. `FALSE` everywhere
else. A marginal likelihood criterion routes on this: where it is `TRUE`
the penalty's contribution to the criterion is a log determinant, which
[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md),
[`penalty_rank()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md),
[`penalty_null_basis()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
and
[`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
supply, and where it is `FALSE` those four reject.

## Usage

``` r
is_quadratic(pen, ...)
```

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object of any branch.

- ...:

  Passed to methods. No shipped method reads it.

## Value

A single logical.

## The quantity behind it

A quadratic penalty carries a fixed matrix \\P\\ and one smoothing
parameter \\\lambda\\, and is the negative log-density of the Gaussian
prior with precision \\\lambda P\\ on \\D\beta\\:

\$\$\rho(\beta; \lambda) = \tfrac{\lambda}{2} (D\beta)^\top P (D\beta) -
\tfrac{1}{2}\log^{+}\lvert \lambda P \rvert + \tfrac{r}{2}\log(2\pi),
\qquad \log^{+}\lvert \lambda P \rvert = r \log \lambda + \log^{+}\lvert
P \rvert,\$\$

with \\\log^{+}\\ the log pseudo-determinant and \\r =
\operatorname{rank}(P)\\. A structured penalty is the same shape with
\\S(\theta)\\ a parameters7 matrix parameter in place of \\\lambda P\\,
so its matrix moves with several hyperparameters and its log-determinant
is not linear in any one of them.

## The three branches it is TRUE for

[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
is the third, its matrix being \\\sum_k \lambda_k P_k\\: a sum of
quadratic forms is a quadratic form, so it has all four quantities. What
is particular about it is that the matrix moves with one hyperparameter
per component rather than with a single scale, which is the property
[`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
routes on when it decides whether to read a slope or to compare a
gradient.

## Methods

The method on the base class
[`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
returns `FALSE`, so a branch that is quadratic says so by registering a
method of its own and every other branch is answered without one.

## See also

[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
and the three quantities beside it,
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
and
[`has_prox()`](https://statmodels7.github.io/penalties7/reference/has_prox.md)
for the other two predicates, and
[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md),
[`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md)
and
[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
for the branches.

## Examples

``` r
# The two branches this is TRUE for.
is_quadratic(quadratic_penalty(diag(2)))
#> [1] TRUE
is_quadratic(ridge_penalty())
#> [1] TRUE
is_quadratic(structured_penalty(
  parameters7::log_cholesky(2, role = "precision")))
#> [1] TRUE

# The third, whose matrix moves with one hyperparameter per component.
add <- additive_penalty(list(diag(3), diag(c(1, 1, 0))))
c(quadratic = is_quadratic(add), rank = penalty_rank(add))
#> quadratic      rank 
#>         1         3 

# A penalty with a kink has none of it.
is_quadratic(lasso_penalty())
#> [1] FALSE
```
