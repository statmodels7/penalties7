# Proximal Operator of a Quadratic or Structured Penalty

One linear solve, \\(I + tS)^{-1}v\\, with \\S\\ the penalty's Hessian:
\\\lambda D'PD\\ for a quadratic penalty and \\\Omega(\theta)\\ for a
structured one. The two methods have the same body, because the
objective is quadratic in both.

## Arguments

- pen:

  A
  [`QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/QuadraticPenalty.md)
  or
  [`StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/StructuredPenalty.md)
  object.

- v:

  A numeric vector of length `pen@n_coef`.

- step:

  The step length \\t\\, a single positive finite number.

- theta:

  A named list of hyperparameter values.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of the same length as `v`.

## Details

The subproblem is \\\arg\min\_\beta \\\tfrac{1}{2t}\lVert\beta -
v\rVert^2 + \tfrac{1}{2}\beta'S\beta\\\\, whose stationary condition is
\\(\beta - v)/t + S\beta = 0\\. The normalizing constant does not depend
on \\\beta\\ and drops out.

The objective stays quadratic under any linear map, so this route places
no restriction on \\D\\; the separable branches need the identity map or
a diagonal one, and a structured penalty has no map. When \\S\\ is a
base matrix the solve is dense and costs \\O(q^3)\\. For a penalty built
with `blocks > 1`, \\S\\ is a `dgCMatrix`, the sum
`diag(length(v)) + step * S` stays sparse, and
[`solve()`](https://rdrr.io/r/base/solve.html) uses the sparse
factorization of Matrix.

## See also

[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
for the generic and the other branches,
[`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
for the matrix solved against.

## Examples

``` r
v <- c(2, 0.3, -1.4)

# With S = lambda I the solve is a plain shrinkage.
penalty_prox(quadratic_penalty(diag(3)), v, 1, list(lambda = 1.5))
#> [1]  0.80  0.12 -0.56
v / (1 + 1.5)
#> [1]  0.80  0.12 -0.56

# A first-difference penalty leaves constant vectors unpenalized, so it
# pulls the coefficients towards their mean and not towards zero.
penalty_prox(quadratic_penalty(crossprod(diff(diag(3)))), v, 1,
             list(lambda = 5))
#> [1] 0.58333333 0.30000000 0.01666667

# A structured penalty at a zero log-Cholesky free vector is the ridge with
# lambda = 1.
s <- structured_penalty(parameters7::log_cholesky(3))
penalty_prox(s, v, 1,
             as.list(stats::setNames(rep(0, 6), s@params)))
#> [1]  1.00  0.15 -0.70
```
