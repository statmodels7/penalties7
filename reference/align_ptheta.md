# Align and Validate the Hyperparameters

Reorders `theta` by name, strips stray names off the values and
validates against `params_bounds` treated as open intervals – the
distributions7 contract, restated here for hyperparameters. A named
numeric vector is accepted in place of the list and converted to one, so
that every branch reads the same shape.

## Usage

``` r
align_ptheta(pen, theta)
```

## Arguments

- pen:

  A
  [`penalty`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- theta:

  A named list of hyperparameter values, or a named numeric vector
  carrying the same.

## Value

The aligned list.
