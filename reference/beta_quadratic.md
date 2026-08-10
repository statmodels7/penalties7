# Is a Penalty Quadratic in the Coefficients?

`TRUE` when \\\partial^3\rho/\partial\beta^3\\ is exactly zero, which
spares a consumer that third derivative altogether.

## Usage

``` r
beta_quadratic(pen, theta, ...)
```

## Arguments

- pen:

  A
  [`penalty`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- theta:

  A named list of hyperparameters.

- ...:

  Passed to methods.

## Value

A single logical.

## Details

Distinct from
[`is_quadratic`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md),
which is about the whole construction of
[`quadratic_penalty`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md):
a structured penalty is quadratic in the coefficients and is not a
quadratic penalty, and a
[`distrib_penalty`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md)
over a gaussian parent is quadratic in the coefficients while its
Hessian is not linear in its hyperparameters. The three properties are
independent and each is asked of the penalty.

## See also

[`penalty_dhessian`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)

## Examples

``` r
beta_quadratic(quadratic_penalty(diag(3)), list(lambda = 1))
#> [1] TRUE
```
