# Check a Penalty Numerically

Compares the derivatives that a penalty declares with numerical
references, and returns one row per comparison carrying the worst scaled
error and a pass or fail. The gradient is compared with a numDeriv
gradient of the value, and the Hessian with a numerical Jacobian of the
analytic gradient. Each hyperparameter gradient is compared with a
numerical derivative of the value in that hyperparameter, each entry of
the hyperparameter Hessian, diagonal or off-diagonal, with a numerical
derivative of the analytic hyperparameter gradient, and each mixed block
with a numerical Jacobian of the analytic coefficient gradient, so that
no reference is a difference of another difference. A quadratic penalty
also has its three-point identity, its log pseudo-determinant, its rank
and its null basis tested. The function serves to verify the derivatives
of a penalty written by the user.

## Usage

``` r
check_penalty(pen, beta = NULL, theta = NULL, tol = 1e-06, verbose = TRUE)
```

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object, of any branch.

- beta:

  A numeric coefficient vector of length `pen@n_coef`. `NULL`, the
  default, draws one from `rnorm(sd = 1.3)` rounded to two places and
  shifted by `0.11`, then moves it away from the kinks as above. The
  draw is taken from a fixed seed, so the report is the same on two
  runs, and the caller's own `.Random.seed` is restored on exit (or
  removed, if there was none), so a call in the middle of a simulation
  leaves that simulation unchanged.

- theta:

  A named list of hyperparameter values, or a named numeric vector
  carrying the same. `NULL`, the default, places each hyperparameter six
  tenths of the way across its own bounds, reading an infinite lower
  bound as `-1` and an infinite upper bound as two above the lower. That
  gives `lambda = 1.2` for a hyperparameter bounded on \\(0, \infty)\\,
  `alpha = 0.6` for the elastic net, `a = 3.2` for SCAD and
  `gamma = 2.2` for MCP.

- tol:

  The scaled error at or above which a row is reported as `FAILED`. A
  single positive number, `1e-6` by default. A row whose error is not
  finite is also reported as `FAILED`, and so is a row whose numerical
  reference cannot be computed.

- verbose:

  `TRUE`, the default, prints the table without row names. The result is
  returned invisibly either way.

## Value

A data frame with one row per check and three columns: `check`
(character, the row's name as tabulated above), `max_error` (numeric)
and `status` (character, `"OK"` or `"FAILED"`). Returned invisibly. For
the rows that compare with a numerical reference, `max_error` is the
worst absolute difference divided by `max(1, max(abs(reference)))`. For
the three-point identity it is the absolute value of the combination
divided by `max(1, abs(rho(beta)))`; for the slope of the log
pseudo-determinant, the absolute difference between the increase and the
rank; for the rank and for the width of the null basis, the absolute
difference between the two counts; for the annihilation by the null
basis, the largest absolute entry of \\S(\theta)N\\ divided by
`max(1, max(abs(S)))`. A row whose reference signals an error, a
numerical derivative or the eigenvalues of the rank row, has `max_error`
equal to `NaN`.

## The rows

Two rows are produced for every penalty, three more for each
hyperparameter, one more for each pair of distinct hyperparameters, and
four or five more when
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
is `TRUE`:

|  |  |
|----|----|
| row | compares |
| `gradient vs numDeriv` | [`penalty_gradient()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md) against a numerical gradient of [`penalty_value()`](https://statmodels7.github.io/penalties7/reference/penalty_value.md) |
| `hessian vs numDeriv on the gradient` | [`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md) against a numerical Jacobian of the analytic gradient, symmetrized |
| `grad_theta[p] vs numDeriv` | [`penalty_grad_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md) against a numerical derivative of the value in hyperparameter `p` |
| `hess_theta[p_p] vs numDeriv` | the diagonal entry of [`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md) against a numerical derivative of the analytic `grad_theta` |
| `hess_theta[p_r] vs numDeriv` | the off-diagonal entry of [`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md) for hyperparameters `p` and `r`, with `p` before `r` in `pen@params`, against a numerical derivative in `r` of the analytic `grad_theta` of `p` |
| `cross[p] vs Richardson on the gradient` | [`penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md) against a numerical Jacobian of the analytic coefficient gradient in `p` |
| `quadratic three-point identity` | \\\rho(2\beta) - 4\rho(\beta) + 3\rho(0)\\, which vanishes when the value is a quadratic form with no linear term |
| `logpdet linear in log lambda with slope r` | raising `lambda` by a factor of \\e\\ raises the value of [`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md) by the rank |
| `logpdet gradient vs numDeriv` | the `grad` element of [`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md) against a numerical derivative of its `value` |
| `rank vs eigenvalues of the matrix` | [`penalty_rank()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md) against the number of eigenvalues of \\S(\theta)\\, or of the sum of the normalized components of an additive penalty, that exceed a relative tolerance times the largest |
| `null basis width vs n_coef minus rank` | the number of columns of [`penalty_null_basis()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md) against `pen@n_coef` minus [`penalty_rank()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md) |
| `null basis annihilates the matrix` | \\S(\theta)N\\, where \\S(\theta)\\ is [`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md) and \\N\\ is [`penalty_null_basis()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md) |

The two log pseudo-determinant rows are alternatives, as described
below. The row of the annihilation appears only when the null basis has
columns. With \\p\\ hyperparameters the count is therefore \\2 + 3p +
p(p-1)/2\\ for a penalty that is not quadratic, and four or five more
for a quadratic one. A lasso gives 5 rows, SCAD 9, a full-rank ridge 9,
a quadratic penalty over second differences 10, and a structured penalty
over a `3 x 3` log-Cholesky precision 39.

## Which log pseudo-determinant row a quadratic branch gets

The quadratic rows are produced for every branch for which
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
returns `TRUE`: the three-point identity, one of the two log
pseudo-determinant rows, the row of the rank, the row of the width of
the null basis, and the row of the annihilation when the null basis has
columns. When `"lambda"` is among `pen@params`, the matrix is one scale
times a constant, the log pseudo-determinant is linear in
\\\log\lambda\\ with slope the rank, and the row tests that identity.
Otherwise the matrix depends on several hyperparameters, as in a
structured or an additive penalty, and the row compares the reported
gradient with a numerical derivative of the value. A two-component
additive penalty whose sum has a null space therefore produces 14 rows:
the 9 rows of its two hyperparameters, the three-point identity, the
gradient of the log pseudo-determinant, the rank, the width of the null
basis and its annihilation. When the sum has full rank the row of the
annihilation is absent and the count is 13.

The row of the rank counts the eigenvalues that exceed a relative
tolerance times the largest. For a quadratic or a structured penalty
they are the eigenvalues of \\S(\theta)\\ at the `theta` of the check.
For an additive penalty they are the eigenvalues of the sum of the
components, each divided by its largest absolute entry, so that the row
reproduces the count of
[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md):
a count on \\S(\theta)\\ drops the smaller components once their scales
differ by more than the tolerance, as that function's page describes.
The tolerance is the `tol` stored in a quadratic or an additive penalty,
the one given to its constructor, and `1e-10` for a structured penalty,
which stores none. The row's `max_error` is the absolute difference
between that count and
[`penalty_rank()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md).
A stored rank that differs from the count by one or more therefore fails
the row of the rank, and a null basis whose number of columns differs
from `pen@n_coef` minus the stored rank fails the row of the width.

## What is not checked

The check covers the value and its derivatives on the parameter scale,
and the quantities that a marginal criterion reads. It does not test the
`hess` element of
[`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md),
the `scale = "link"` option of the hyperparameter derivatives,
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md),
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md),
[`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md),
[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md),
[`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md),
[`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md),
[`penalty_dhessian_beta()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian_beta.md),
[`penalty_d2hessian_beta()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.md),
[`penalty_dhessian_beta_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.md),
[`penalty_draw()`](https://statmodels7.github.io/penalties7/reference/penalty_draw.md)
or
[`penalty_readable()`](https://statmodels7.github.io/penalties7/reference/penalty_readable.md).
It calls
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
to place the coefficients but does not test the kinks themselves.

## Where the coefficients are placed

A penalty with a kink has no derivative there, so a numerical reference
that straddles one measures the kink and not the formula. With `beta`
left at `NULL`, `0.033` is added to every coefficient, up to fifty
times, until every coordinate of \\D\beta\\ lies more than `0.05` from
every kink that
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
returns. A draw that already clears the kinks is used unchanged, and
after fifty steps the last draw is used whether or not it clears them.

## Errors

numDeriv is in `Suggests` and the derivative rows need it, so the
function signals an error when it is not installed. A `theta` outside
the penalty's open bounds or missing a hyperparameter is rejected before
any check runs. The analytic methods are evaluated at `beta` and `theta`
outside the numerical references, so an error in
[`penalty_value()`](https://statmodels7.github.io/penalties7/reference/penalty_value.md),
[`penalty_gradient()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md),
[`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md),
the hyperparameter derivatives or the quadratic quantities is signaled
as it arises. An error inside a numerical reference is not signaled:
numDeriv stops when the function it differentiates returns a value that
is not finite, and the row whose reference stopped is then reported with
`max_error` equal to `NaN` and status `FAILED`, while the other rows are
computed as usual.

## See also

[`penalty_value()`](https://statmodels7.github.io/penalties7/reference/penalty_value.md)
and
[`penalty_gradient()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
for the quantities checked,
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
for the set that the coefficients avoid,
[`numericals7::check_abs_smoother()`](https://statmodels7.github.io/numericals7/reference/check_abs_smoother.html)
for the same check on a smoother,
[`linkfunctions7::check_link()`](https://statmodels7.github.io/linkfunctions7/reference/check_link.html)
and
[`distributions7::check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.html)
for the corresponding checks of a link and of a distribution.

## Examples

``` r
# First differences over five coefficients: rank 4, one null direction.
pen <- quadratic_penalty(crossprod(diff(diag(5))))
res <- check_penalty(pen)
#>                                        check    max_error status
#>                         gradient vs numDeriv 4.007727e-11     OK
#>          hessian vs numDeriv on the gradient 7.928288e-12     OK
#>               grad_theta[lambda] vs numDeriv 1.851376e-11     OK
#>        hess_theta[lambda_lambda] vs numDeriv 1.488393e-11     OK
#>  cross[lambda] vs Richardson on the gradient 1.291110e-11     OK
#>               quadratic three-point identity 2.296653e-16     OK
#>    logpdet linear in log lambda with slope r 4.440892e-16     OK
#>            rank vs eigenvalues of the matrix 0.000000e+00     OK
#>        null basis width vs n_coef minus rank 0.000000e+00     OK
#>            null basis annihilates the matrix 2.775558e-16     OK
all(res$status == "OK")
#> [1] TRUE

# A separable penalty has fewer rows: no matrix, so no quadratic checks.
nrow(check_penalty(lasso_penalty(n_coef = 4), verbose = FALSE))
#> [1] 5

# A penalty with a wrong gradient fails the rows that read the gradient.
# The broken method is registered on a subclass, because registering it on
# the real class would change the generic for the rest of the session.
Broken <- S7::new_class("Broken", parent = QuadraticPenalty)
bad <- do.call(Broken, S7::props(quadratic_penalty(diag(3))))
S7::method(penalty_gradient, Broken) <- function(pen, beta, theta, ...) {
  1.05 * S7::method(penalty_gradient, QuadraticPenalty)(pen, beta, theta, ...)
}
failed <- check_penalty(bad, verbose = FALSE)
failed[failed$status != "OK", c("check", "max_error")]
#>                                         check  max_error
#> 1                        gradient vs numDeriv 0.05000000
#> 2         hessian vs numDeriv on the gradient 0.04761905
#> 5 cross[lambda] vs Richardson on the gradient 0.04761905
```
