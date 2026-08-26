# A Zero Entry for Every Hyperparameter Pair

Returns a list holding the same zero object under every
hyperparameter-pair key, which is the answer of a penalty whose Hessian
is linear in its hyperparameters.

## Usage

``` r
zero_pairs(pen, z)
```

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- z:

  The zero object: a `pen@n_coef` by `pen@n_coef` matrix for
  [`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md),
  a numeric vector of that length for
  [`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md).

## Value

A named list of \\p(p+1)/2\\ entries, keyed as
[`ptheta_pairs()`](https://statmodels7.github.io/penalties7/reference/ptheta_pairs.md)
keys them, each holding `z`.

## Details

The quadratic and additive branches use it for both
[`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md)
and
[`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md),
the two differing only in the shape of the zero.

## See also

[`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md),
[`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md),
[`ptheta_pairs()`](https://statmodels7.github.io/penalties7/reference/ptheta_pairs.md)
