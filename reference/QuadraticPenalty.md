# S7 Class for the Quadratic Penalty

The class
[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
builds. Beyond the seven properties every penalty carries it stores the
matrix \\P\\, its rank, an orthonormal basis of the null space of
\\D'PD\\, the log pseudo-determinant of \\P\\ and the assembled
\\D'PD\\. All five are fixed by one eigendecomposition at construction
and none of them moves with the hyperparameter, so every quantity the
branch supplies is closed form in \\\lambda\\.

## Usage

``` r
QuadraticPenalty(
  penalty_name = character(0),
  map = NULL,
  n_coef = integer(0),
  params = character(0),
  params_bounds = list(),
  link_params = list(),
  params_smooth = logical(0),
  P = NULL,
  p_rank = integer(0),
  null_basis = NULL,
  logpdet_P = integer(0),
  DPD = NULL
)
```

## Arguments

- penalty_name:

  A single string naming the penalty, used by
  [`print()`](https://rdrr.io/r/base/print.html) and by consumers that
  report which penalty a block carries.

- map:

  The matrix \\D\\, of \\m\\ rows and \\q\\ columns, or `NULL` for the
  identity. A Matrix object is kept in its own storage; a diagonal map
  is what standardization comes to and is recognized by its class.

- n_coef:

  The number of coefficients \\q\\. A single whole number.

- params:

  The hyperparameter names, in the order every derivative list is keyed
  by. `character(0)` for a penalty with none.

- params_bounds:

  A named list, one entry per hyperparameter, each a numeric pair giving
  an **open** interval. A value at either endpoint is rejected, so
  `(0, Inf)` excludes zero.

- link_params:

  A named list, one linkfunctions7 link per hyperparameter, carrying
  that hyperparameter's own interval onto the whole real line.

- params_smooth:

  A logical vector, one entry per hyperparameter, `TRUE` where the value
  is differentiable in it.

- P:

  The symmetric positive semidefinite matrix of the quadratic form, with
  as many rows as the map has, or as many as there are coefficients when
  the map is `NULL`. A `dgCMatrix` under `blocks > 1`.

- p_rank:

  The rank of \\P\\, a single whole number, counted at construction by a
  relative eigenvalue rule.

- null_basis:

  An orthonormal basis of the null space of \\D'PD\\, with `n_coef` rows
  and `n_coef - p_rank` columns, and no columns when the penalty is full
  rank.

- logpdet_P:

  The log pseudo-determinant of \\P\\: the sum of the logarithms of its
  non-zero eigenvalues. A single number.

- DPD:

  The assembled \\D'PD\\, an `n_coef` by `n_coef` matrix, cached so that
  the gradient and the Hessian never re-form it.

## Value

An S7 object of class `QuadraticPenalty`, inheriting from
[`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md),
with the seven inherited properties and the five above.

## Details

Construct through
[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md),
which computes the five derived properties and checks that \\P\\ is
symmetric. This raw constructor takes them as given and checks nothing,
so a wrong rank here would be reported by
[`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
and by nothing else.

## See also

[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
for the constructor to use,
[`penalty_value.QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.QuadraticPenalty.md)
for what the branch computes,
[`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md)
for the other branch that answers
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
with `TRUE`.

## Examples

``` r
pen <- quadratic_penalty(crossprod(diff(diag(5), differences = 2)))
S7::S7_inherits(pen, QuadraticPenalty)
#> [1] TRUE

# The five derived properties, fixed at construction.
pen@p_rank
#> [1] 3
dim(pen@null_basis)
#> [1] 5 2
pen@logpdet_P
#> [1] 3.912023
```
