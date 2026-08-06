# Value of a Penalty

The scalar \\\rho(D\beta;\theta)\\, with the normalizing constant
included whenever the penalty is proper, so that the value is exactly
the negative log-density of the prior.

## Usage

``` r
penalty_value(pen, beta, theta, ...)
```

## Arguments

- pen:

  A
  [`penalty`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- beta:

  A numeric vector of coefficients.

- theta:

  A named list of hyperparameter values.

- ...:

  Passed to methods.

## Value

A single number.

## Examples

``` r
pen <- quadratic_penalty(diag(3))
penalty_value(pen, c(1, 0, -1), list(lambda = 2))
#> [1] 3.717095
```
