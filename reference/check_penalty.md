# Check a Penalty Numerically

Compares every closed form a penalty declares against a route that
shares no code with it, and returns one row per comparison carrying the
worst relative error and a pass or fail. The gradient and the Hessian go
against numDeriv on the value; each hyperparameter block goes against
numDeriv in that hyperparameter, the mixed block by Richardson on the
analytic gradient so that no difference is ever taken of another
difference; and a quadratic penalty has its three-point identity, its
log pseudo-determinant and its null basis tested as well. Write a
penalty of your own and this says whether its derivatives are right.

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
  shifted by `0.11`, then pushes it clear of the kinks as above. The
  draw is taken from a fixed seed, so the report is the same on two
  runs, and the caller's own `.Random.seed` is restored on exit, so a
  call in the middle of a simulation leaves that simulation unchanged.

- theta:

  A named list of hyperparameter values, or a named numeric vector
  carrying the same. `NULL`, the default, places each hyperparameter six
  tenths of the way across its own bounds, reading an infinite lower
  bound as `-1` and an infinite upper bound as two above the lower. That
  gives `lambda = 1.2` on every branch that has one, `alpha = 0.6` for
  the elastic net, `a = 3.2` for SCAD and `gamma = 2.2` for MCP.

- tol:

  The relative error above which a row is reported as `FAILED`. A single
  positive number, `1e-6` by default. Over the nine shipped branches the
  worst error measured is `2.1e-10`, four orders under it, so the
  default separates a correct penalty from one whose formula is wrong in
  its fifth digit.

- verbose:

  `TRUE`, the default, prints the table without row names. The result is
  returned invisibly either way.

## Value

A data frame with one row per check and three columns: `check`
(character, the row's name as tabulated above), `max_error` (numeric,
the worst absolute difference divided by `max(1, max(abs(reference)))`),
and `status` (character, `"OK"` or `"FAILED"`). Returned invisibly.

## The rows

Two rows are produced for every penalty, three more for each
hyperparameter, and up to three more when
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
is `TRUE`:

|  |  |
|----|----|
| row | compares |
| `gradient vs numDeriv` | [`penalty_gradient()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md) against a numerical gradient of [`penalty_value()`](https://statmodels7.github.io/penalties7/reference/penalty_value.md) |
| `hessian vs numDeriv on the gradient` | [`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md) against a numerical Jacobian of the analytic gradient, symmetrized |
| `grad_theta[p] vs numDeriv` | [`penalty_grad_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md) against a numerical derivative of the value in hyperparameter `p` |
| `hess_theta[p_p] vs numDeriv` | [`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md) against a numerical derivative of the analytic `grad_theta` |
| `cross[p] vs Richardson on the gradient` | [`penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md) against a numerical Jacobian of the analytic coefficient gradient in `p` |
| `quadratic three-point identity` | \\\rho(2\beta) - 4\rho(\beta) + 3\rho(0)\\, which vanishes when the value is a quadratic form with no linear term |
| `logpdet linear in log lambda with slope r` | raising `lambda` by a factor of \\e\\ raises [`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md) by exactly the rank |
| `logpdet gradient vs numDeriv` | the `grad` element of [`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md) against a numerical derivative of its `value`, for a quadratic penalty whose hyperparameters are not a single `lambda` |
| `null basis annihilates the matrix` | \\PN\\, where `N` is [`penalty_null_basis()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md) |

The last two of the log-determinant rows are alternatives: the first is
taken when `lambda` is among the hyperparameters, which is the plain
quadratic branch of one scale multiplying a constant matrix, and the
second otherwise, which is the structured branch. The null-basis row
appears only when the null basis has columns.

So the count is `2 + 3 * length(pen@params)` plus one, two or three. A
lasso gives 5 rows, a full-rank ridge 7, a quadratic penalty over second
differences 8, and a structured penalty over a `3 x 3` log-Cholesky
precision 22.

## What is not checked

The pass is over the value and its derivatives, and over the pieces a
marginal criterion reads. It does not touch
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md),
[`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md),
[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md),
[`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md),
[`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md)
or
[`penalty_readable()`](https://statmodels7.github.io/penalties7/reference/penalty_readable.md).
It calls
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
to place the grid but never tests the kinks themselves.

## Which logpdet row a quadratic branch gets

The three quadratic rows run for every branch
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
answers `TRUE` for, and the middle one comes in two shapes. Where the
matrix is one scale times a constant, which is
`"lambda" %in% pen@params`, the log pseudo-determinant is linear in
`log lambda` with slope the rank, and that identity is read directly.
Where it moves with several hyperparameters, as in a structured or an
additive penalty, there is no such slope and the row compares the
reported gradient against `numDeriv` instead. A two-component additive
penalty therefore produces 11 rows: the 8 its two hyperparameters earn,
plus the three-point identity, the log pseudo-determinant's gradient and
the null basis.

## Where the grid is placed

A penalty with a kink has no derivative there, so a numerical reference
straddling one measures the kink and not the formula. With `beta` left
at `NULL` the draw is nudged in steps of `0.033`, up to fifty times,
until every coordinate of \\D\beta\\ sits at least `0.05` from every
kink the object declares. At the default coefficient count the first
draw already clears the SCAD and MCP kink sets and no step is taken.

## Errors

numDeriv is in `Suggests` and every row needs it, so the function stops
when it is not installed. A `theta` outside the penalty's open bounds or
missing a hyperparameter is rejected before any check runs.

## See also

[`penalty_value()`](https://statmodels7.github.io/penalties7/reference/penalty_value.md)
and
[`penalty_gradient()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
for the quantities checked,
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
for the set the grid avoids,
[`check_abs_smoother()`](https://statmodels7.github.io/penalties7/reference/check_abs_smoother.md)
for the same service on a smoother,
[`linkfunctions7::check_link()`](https://statmodels7.github.io/linkfunctions7/reference/check_link.html)
and
[`distributions7::check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.html)
for the siblings this follows.

## Examples

``` r
# Second differences over five coefficients: rank 4, one null direction.
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
#>            null basis annihilates the matrix 2.775558e-16     OK
all(res$status == "OK")
#> [1] TRUE

# A separable penalty has fewer rows: no matrix, so no quadratic checks.
nrow(check_penalty(lasso_penalty(n_coef = 4), verbose = FALSE))
#> [1] 5

# The validator earns its keep on a penalty that is wrong. Register the
# broken method on a subclass: registering on the real class would mutate
# the generic for the rest of the session.
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
