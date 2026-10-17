# S7 Class for the Quadratic Penalty

`QuadraticPenalty` is the class of the objects that
[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
builds. Beyond the seven properties that every penalty carries, it
stores the matrix \\P\\, the rank of \\D'PD\\, an orthonormal basis of
the null space of \\D'PD\\, the log pseudo-determinant of \\D'PD\\ and
the assembled \\D'PD\\. All five are computed once at construction and
are independent of the hyperparameter, so every quantity of the branch
is closed form in \\\lambda\\.

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
  logpdet_DPD = integer(0),
  DPD = NULL,
  tol = 1e-10
)
```

## Arguments

- penalty_name:

  A single string naming the penalty, used by
  [`print()`](https://rdrr.io/r/base/print.html) and by any code that
  reports which penalty a block carries.

- map:

  The matrix \\D\\, of \\m\\ rows and \\q\\ columns, or `NULL` for the
  identity. A Matrix object is kept in its own storage, so a diagonal
  map, such as the map that standardization produces, is recognized as
  diagonal by its class.

- n_coef:

  The number of coefficients \\q\\. The value is a single whole number.

- params:

  The hyperparameter names, in the order by which every derivative list
  is keyed. The value is `character(0)` for a penalty with none.

- params_bounds:

  A named list with one entry per hyperparameter, and each entry is a
  numeric pair giving an open interval. A value at either endpoint is
  rejected, so `(0, Inf)` excludes zero.

- link_params:

  A named list, one linkfunctions7 link per hyperparameter, carrying
  that hyperparameter's own interval onto the whole real line.

- params_smooth:

  A logical vector, one entry per hyperparameter, `TRUE` where the value
  is differentiable in it.

- P:

  The symmetric positive semidefinite matrix of the quadratic form, with
  as many rows as the map has, or as many as there are coefficients when
  the map is `NULL`. It is a `dgCMatrix` under `blocks > 1`.

- p_rank:

  The rank of \\D'PD\\, a single whole number, counted at construction
  by a relative eigenvalue rule. Without a map it is the rank of \\P\\.

- null_basis:

  An orthonormal basis of the null space of \\D'PD\\, with `n_coef` rows
  and `n_coef - p_rank` columns, and no columns when the penalty is full
  rank.

- logpdet_DPD:

  The log pseudo-determinant of \\D'PD\\: the sum of the logarithms of
  its non-zero eigenvalues. It is a single number. Without a map it is
  the log pseudo-determinant of \\P\\.

- DPD:

  The assembled \\D'PD\\, an `n_coef` by `n_coef` matrix, cached so that
  the gradient and the Hessian never re-form it.

- tol:

  The relative eigenvalue tolerance, a single positive number, `1e-10`
  by default.
  [`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
  counts `p_rank` at this tolerance and stores it, and
  [`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
  counts the eigenvalues of the matrix at the same tolerance when it
  compares them with the rank.

## Value

An S7 object of class `QuadraticPenalty`, inheriting from
[`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md),
with the seven inherited properties and the six above.

## Details

Construct through
[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md),
which computes the five properties and checks that \\P\\ is symmetric
and positive semidefinite. This raw constructor takes them as given and
does not check any of them. A wrong rank passed here enters the value,
the hyperparameter derivatives,
[`penalty_rank()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md),
[`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
and
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md),
and
[`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
does not detect it, because every one of its rows reads the same stored
rank.

## Methods

This section lists the methods registered on this class, by the page
that documents each.

- The value and its coefficient derivatives:
  [penalty_value()](https://statmodels7.github.io/penalties7/reference/penalty_value.QuadraticPenalty.md),
  [penalty_gradient()](https://statmodels7.github.io/penalties7/reference/penalty_gradient.QuadraticPenalty.md)
  and
  [penalty_hessian()](https://statmodels7.github.io/penalties7/reference/penalty_gradient.QuadraticPenalty.md).

- The hyperparameter derivatives:
  [penalty_grad_theta()](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.QuadraticPenalty.md),
  [penalty_hess_theta()](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.QuadraticPenalty.md)
  and
  [penalty_cross()](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.QuadraticPenalty.md).

- The predicates:
  [penalty_kinks()](https://statmodels7.github.io/penalties7/reference/penalty_kinks.QuadraticPenalty.md),
  [is_proper()](https://statmodels7.github.io/penalties7/reference/penalty_kinks.QuadraticPenalty.md)
  and
  [is_quadratic()](https://statmodels7.github.io/penalties7/reference/penalty_kinks.QuadraticPenalty.md).

- The quantities of a marginal criterion:
  [penalty_matrix()](https://statmodels7.github.io/penalties7/reference/penalty_matrix.QuadraticPenalty.md),
  [penalty_rank()](https://statmodels7.github.io/penalties7/reference/penalty_matrix.QuadraticPenalty.md),
  [penalty_null_basis()](https://statmodels7.github.io/penalties7/reference/penalty_matrix.QuadraticPenalty.md)
  and
  [penalty_logpdet()](https://statmodels7.github.io/penalties7/reference/penalty_matrix.QuadraticPenalty.md).

- The higher derivatives of a marginal criterion:
  [penalty_dhessian()](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.QuadraticPenalty.md),
  [penalty_d2hessian()](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.QuadraticPenalty.md),
  [penalty_dcross()](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.QuadraticPenalty.md)
  and
  [beta_quadratic()](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.QuadraticPenalty.md).

- The movement of the Hessian with the coefficients:
  [penalty_dhessian_beta()](https://statmodels7.github.io/penalties7/reference/penalty_dhessian_beta.penalty.md),
  [penalty_d2hessian_beta()](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.penalty.md)
  and
  [penalty_dhessian_beta_theta()](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.penalty.md).
  The last two return zero matrices because the Hessian does not depend
  on the coefficients, and are documented on the page of the base-class
  methods.

- The proximal operator:
  [penalty_prox()](https://statmodels7.github.io/penalties7/reference/penalty_prox.QuadraticPenalty.md).

The class has no
[`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md),
[`penalty_draw()`](https://statmodels7.github.io/penalties7/reference/penalty_draw.md)
or
[`penalty_readable()`](https://statmodels7.github.io/penalties7/reference/penalty_readable.md)
method of its own, so those generics dispatch to the methods of the base
class
[`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md).

## See also

[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
for the constructor to use,
[`penalty_value.QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.QuadraticPenalty.md)
for what the branch computes,
[`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md)
and
[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
for the other branches for which
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
returns `TRUE`.

## Examples

``` r
pen <- quadratic_penalty(crossprod(diff(diag(5), differences = 2)))
S7::S7_inherits(pen, QuadraticPenalty)
#> [1] TRUE

# These are three of the six properties, all computed at construction.
pen@p_rank
#> [1] 3
dim(pen@null_basis)
#> [1] 5 2
pen@logpdet_DPD
#> [1] 3.912023
```
