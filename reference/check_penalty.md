# Check a Penalty Numerically

The sibling of `check_link` and `check_distrib`: every closed form is
compared against a route that shares no code with it. The gradient and
the Hessian are checked against numDeriv on the value, the theta blocks
against numDeriv in each hyperparameter – the mixed block by Richardson
on the analytic gradient, never a nested difference – and the map by
comparing \\\rho(D\beta)\\ routes. Grids are placed away from the kink
set the object itself declares.

## Usage

``` r
check_penalty(pen, beta = NULL, theta = NULL, tol = 1e-06, verbose = TRUE)
```

## Arguments

- pen:

  A
  [`penalty`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- beta:

  A coefficient vector, or `NULL` for a default draw.

- theta:

  A named hyperparameter list, or `NULL` for midpoints.

- tol:

  The comparison tolerance.

- verbose:

  Logical; print the table.

## Value

A data frame with one row per check, invisibly when printed.

## See also

[`penalty_value`](https://statmodels7.github.io/penalties7/reference/penalty_value.md),
[`has_prox`](https://statmodels7.github.io/penalties7/reference/has_prox.md)

## Examples

``` r
res <- check_penalty(quadratic_penalty(diag(3)))
#>                                        check    max_error status
#>                         gradient vs numDeriv 1.860032e-11     OK
#>          hessian vs numDeriv on the gradient 1.076602e-11     OK
#>               grad_theta[lambda] vs numDeriv 5.633756e-12     OK
#>        hess_theta[lambda_lambda] vs numDeriv 3.183821e-11     OK
#>  cross[lambda] vs Richardson on the gradient 3.355036e-12     OK
#>               quadratic three-point identity 1.810552e-16     OK
#>    logpdet linear in log lambda with slope r 4.440892e-16     OK
all(res$status == "OK")
#> [1] TRUE
```
