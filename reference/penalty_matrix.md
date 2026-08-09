# The Pieces a Marginal Criterion Consumes

For a quadratic penalty: `penalty_matrix` returns \\\lambda D'PD\\,
`penalty_rank` its rank (fixed at construction), `penalty_null_basis`
the exact null basis for the model layer to intersect across terms, and
`penalty_logpdet` the log pseudo-determinant \\r\log\lambda +
\log\mathrm{pdet}(P)\\ with its first two theta derivatives. Every other
penalty rejects: a marginal criterion for a non-Gaussian prior is not a
determinant, and pretending otherwise would produce numbers silently.

## Usage

``` r
penalty_matrix(pen, theta, ...)

penalty_rank(pen, ...)

penalty_null_basis(pen, ...)

penalty_logpdet(pen, theta, ...)
```

## Arguments

- pen:

  A
  [`penalty`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- theta:

  A named list of hyperparameter values.

- ...:

  Passed to methods.

## Value

`penalty_matrix` a `q x q` matrix; `penalty_rank` an integer;
`penalty_null_basis` a matrix with `q` rows (zero columns when the
penalty is full rank); `penalty_logpdet` a list with elements `value`,
`grad` and `hess`.

## See also

[`is_quadratic`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md),
[`quadratic_penalty`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md),
[`additive_penalty`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)

## Examples

``` r
pen <- quadratic_penalty(crossprod(diff(diag(4))))
penalty_rank(pen)
#> [1] 3
penalty_logpdet(pen, list(lambda = 2))$value
#> [1] 3.465736
```
