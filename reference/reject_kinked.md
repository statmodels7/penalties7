# Reject a Kinked Penalty

Signals an error naming the penalty and the generic when the penalty
carries a kink, and returns invisibly otherwise. Called at the head of
each of the separable branch's three marginal derivative methods.

## Usage

``` r
reject_kinked(pen, what)
```

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- what:

  The generic's name, a single string, quoted back in the message.

## Value

`NULL`, invisibly, when there is no kink.

## Details

A marginal criterion is a Laplace expansion at the penalized mode, and a
kinked penalty puts coefficients exactly at the kink, where the third
derivative does not exist. Rejecting there is better than returning the
one-sided value, which would give a criterion no one could interpret.

The kinks are read straight off `pen@kinks` rather than through
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md),
which takes a `theta` this is asked before there is one to pass. An
object without that property, which no shipped branch is, is treated as
having none.

## Errors

`'<penalty>' has a kink, so <what>() does not exist there and its hyperparameters cannot be estimated by a marginal criterion.`

## See also

[`penalty_dhessian.DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.DistribPenalty.md),
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
