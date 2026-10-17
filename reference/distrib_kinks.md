# Where the Parent's Log-Density Has a Kink

Returns the points at which a separable penalty built on this
distribution is not differentiable, derived from the distribution
itself, or `numeric(0)` when there are none. Among the families that
distributions7 ships, only the two Laplace parametrizations and the
elastic net have a location in which the log-density is not smooth, so
only these give a kink, and only when
[`distributions7::fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.html)
holds the location.

## Usage

``` r
distrib_kinks(d)
```

## Arguments

- d:

  A univariate distributions7 object. An object without
  `parent_distrib`, or whose chain of parents has no level carrying
  `fixed_params`, yields `numeric(0)` without further inspection.

## Value

A numeric vector of the confirmed kink points, without duplicates and
possibly empty.

## How a candidate is found

A distribution records in `params_smooth` the parameters in which its
log-likelihood is differentiable, and for a location family a location
that is not smooth is a kink in the argument at that location. A penalty
is the negative log-density read at the coefficient, so a
[`distributions7::fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.html)
wrapper that holds such a parameter at a value puts the kink at that
value: `fixed(laplace_distrib(), mu = 0)` is a Laplace prior at zero and
has a kink at zero, while `fixed(laplace_distrib(), mu = 2)` has one at
two.

A free parameter gives no candidate, because its value is set by the
hyperparameters at the time of each evaluation. A Laplace with a free
location therefore has no reported kink.

The wrapper is recognized by the properties that it carries,
`parent_distrib` and `fixed_params`, and not by its class, because the
classes of the
[`distributions7::fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.html)
wrapper are not exported and a family written outside distributions7 may
hold its parameters in the same way. The chain of `parent_distrib` is
walked from the outermost object down, and every level that carries
`fixed_params` contributes candidates, so that
`truncated(fixed(laplace2_distrib(), mu = 0), lower = -1, upper = 3)`
has its kink at zero.

## How a candidate is confirmed

A candidate is dropped when it is not finite or when it does not lie
strictly inside the support of the outermost object, because a point at
the edge of the support is not a kink of the density on its support: a
`folded()` Laplace held at zero has no kink at zero. Each remaining
candidate is then tested by comparing the one-sided derivatives of the
log-density across it at a probe point
([`probe_theta()`](https://statmodels7.github.io/penalties7/reference/probe_theta.md),
[`has_jump()`](https://statmodels7.github.io/penalties7/reference/has_jump.md)),
and a candidate at which the derivative does not jump is dropped. The
inference from `params_smooth` supplies the points to test, and the test
removes the candidates at which the density is smooth, for example a
candidate that comes from a non-smooth parameter that is not a location.
The test is made on the outermost object, so a wrapper that moves the
point drops the candidate.

## See also

[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md),
which calls this function when its argument `kinks` is `NULL`,
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
for the kinks of a built penalty,
[`has_jump()`](https://statmodels7.github.io/penalties7/reference/has_jump.md)
for the test.

## Examples

``` r
# A Laplace held at a location has a kink at that location.
distrib_kinks(distributions7::fixed(distributions7::laplace_distrib(),
                                    mu = 0))
#> [1] 0
distrib_kinks(distributions7::fixed(distributions7::laplace_distrib(),
                                    mu = 2))
#> [1] 2

# A Gaussian is smooth everywhere, so there is no candidate.
distrib_kinks(distributions7::fixed(distributions7::gaussian1_distrib(),
                                    mu = 0))
#> numeric(0)

# A free location gives no candidate.
distrib_kinks(distributions7::laplace_distrib())
#> numeric(0)

# A fixed() level under another wrapper is found; the origin of a folded
# Laplace is the edge of its support and is dropped.
lap0 <- distributions7::fixed(distributions7::laplace2_distrib(), mu = 0)
distrib_kinks(distributions7::truncated(lap0, lower = -1, upper = 3))
#> [1] 0
distrib_kinks(distributions7::folded(lap0))
#> numeric(0)
```
