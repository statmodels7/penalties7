# The Structure's Free Vector From the Aligned Hyperparameters

Unlists the aligned hyperparameter list into the numeric vector the
structure's own generics take, in `pen@params` order, which is the
structure's `free_names` order.

## Usage

``` r
struct_eta(pen, theta)
```

## Arguments

- pen:

  A
  [`StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/StructuredPenalty.md)
  object.

- theta:

  The aligned hyperparameter list, as
  [`align_ptheta()`](https://statmodels7.github.io/penalties7/reference/align_ptheta.md)
  returns it.

## Value

An unnamed numeric vector of length `length(pen@params)`.

## See also

[`struct_omega()`](https://statmodels7.github.io/penalties7/reference/struct_omega.md),
[`align_ptheta()`](https://statmodels7.github.io/penalties7/reference/align_ptheta.md)
