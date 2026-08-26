# Where the Parent's Log-Density Has a Kink

Returns the points at which a separable penalty built on this
distribution would not be differentiable, derived from the distribution
itself. Returns `numeric(0)` where there are none, which is the answer
for every family the toolkit ships except the Laplace and the elastic
net.

## Usage

``` r
distrib_kinks(d)
```

## Arguments

- d:

  A univariate distributions7 object. An object carrying neither
  `parent_distrib` nor `fixed_params` gets `numeric(0)` without further
  inspection.

## Value

A numeric vector of the confirmed kink points, possibly empty, with
duplicates and non-finite candidates removed.

## How a candidate is found

A distribution records which of its parameters the log-likelihood is
differentiable in, through `params_smooth`, and for a location family a
location that is not smooth is a kink in the argument at that location.
A penalty is the negative log-density read in the coefficient, so a
[`distributions7::fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.html)
wrapper holding such a parameter at a value puts the kink there:
`fixed(laplace_distrib(), mu = 0)` is the lasso and has a kink at zero,
while `fixed(laplace_distrib(), mu = 2)` has one at two.

Nothing is taken from a parameter that is **free**, its value being
whatever the hyperparameters say at the time. A Laplace with an unfixed
location therefore reports no kink at all.

## And how it is confirmed

Each candidate is the value its parameter is held at, and whether it is
a kink is then measured rather than inferred, by comparing the one-sided
derivatives of the log-density across it at a probe point
([`probe_theta()`](https://statmodels7.github.io/penalties7/reference/probe_theta.md),
[`has_jump()`](https://statmodels7.github.io/penalties7/reference/has_jump.md)).
Inferring alone would put a kink on any family whose non-smooth
parameter is not a location; measuring alone would need somewhere to
look. A candidate whose derivative does not jump is dropped.

The wrapper is recognized by the properties it carries, `parent_distrib`
and `fixed_params`, and not by its class: the two `fixed` classes are
not exported, and a family written outside distributions7 may hold its
parameters the same way.

## See also

[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md),
which calls this unless told the kinks directly,
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
for what a built penalty reports,
[`has_jump()`](https://statmodels7.github.io/penalties7/reference/has_jump.md)
for the measurement.

## Examples

``` r
# The lasso's parent: a kink at the location it is held at.
distrib_kinks(distributions7::fixed(distributions7::laplace_distrib(),
                                    mu = 0))
#> [1] 0
distrib_kinks(distributions7::fixed(distributions7::laplace_distrib(),
                                    mu = 2))
#> [1] 2

# A Gaussian is smooth everywhere, so there is no candidate to confirm.
distrib_kinks(distributions7::fixed(distributions7::gaussian1_distrib(),
                                    mu = 0))
#> numeric(0)

# And a free location gives nothing, its value not being settled.
distrib_kinks(distributions7::laplace_distrib())
#> numeric(0)
```
