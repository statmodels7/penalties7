# Carry Theta Derivatives Onto the Link Scale

The order 1-2 chain rule with the diagonal Jacobian of the links: the
same interception distributions7 applies, restricted to the two orders a
penalty consumer needs.

## Usage

``` r
ptheta_to_link(pen, theta, g = NULL, H = NULL, cross = NULL)
```

## Arguments

- pen:

  A
  [`penalty`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- theta:

  The aligned hyperparameters.

- g:

  The parameter-scale gradient list, or `NULL`.

- H:

  The parameter-scale Hessian list, or `NULL`.

- cross:

  The parameter-scale mixed list, or `NULL`.

## Value

Whichever of the three was supplied, transformed.
