# A Second Hyperparameter Derivative of a Response Quantity

Differences a parent's analytic mixed derivative once more in each
hyperparameter and carries the result through the penalty's map.

## Usage

``` r
ptheta_second(pen, theta, inner, carry)
```

## Arguments

- pen:

  A
  [`penalty`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- theta:

  A named list of hyperparameters.

- inner:

  A function of `theta` returning a named list, one entry per
  hyperparameter.

- carry:

  A function placing one such entry into coefficient space.

## Value

A named list keyed by hyperparameter pair.

## Details

One central difference of an analytic quantity, which is the fallback
distributions7 uses throughout; the two derivatives act on different
hyperparameters wherever the pair is off the diagonal, and on the same
one on it, where the reference is still analytic and only one difference
is taken.
