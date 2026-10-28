# Reject a Kinked Penalty

Signals an error naming the penalty and the generic when `penalty@kinks`
holds at least one finite value, and returns invisibly otherwise. It is
called at the head of the six methods of the separable branch that need
a third or higher derivative:
[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md),
[`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md),
[`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md),
[`penalty_dhessian_beta()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian_beta.md),
[`penalty_d2hessian_beta()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.md)
and
[`penalty_dhessian_beta_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.md).

## Usage

``` r
reject_kinked(penalty, what)
```

## Arguments

- penalty:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- what:

  The name of the generic, a single string, inserted in the message
  followed by `()`.

## Value

`NULL`, invisibly, when there is no kink.

## Details

A marginal criterion is a Laplace expansion at the penalized mode, and a
kinked penalty places coefficients exactly at the kink, where the third
derivative does not exist. The function therefore signals an error and
does not return a one-sided derivative.

The kinks are read directly from `penalty@kinks` and not through
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md),
because that generic takes a `theta` and this function receives only the
penalty and the name of the generic. A penalty without a `kinks`
property is treated as having no kink; every caller in the package is a
method of
[`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md),
which has the property.

## Errors

`'<penalty>' has a kink, so <what>() does not exist there and its hyperparameters cannot be estimated by a marginal criterion.`

## See also

[`penalty_dhessian.DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.DistribPenalty.md),
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
