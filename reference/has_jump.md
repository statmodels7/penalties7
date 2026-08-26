# Does the Log-Density's Slope Jump Across a Point?

Compares the derivative of the log-density in its argument just either
side of a candidate point, and answers `TRUE` when the two differ by
more than a relative threshold. This is the measurement that turns a
candidate kink, inferred from `params_smooth`, into a confirmed one.

## Usage

``` r
has_jump(d, theta, at, eps = 1e-05)
```

## Arguments

- d:

  A distributions7 object.

- theta:

  A named list of parameter values to read it at, as
  [`probe_theta()`](https://statmodels7.github.io/penalties7/reference/probe_theta.md)
  returns.

- at:

  The candidate point, a single number.

- eps:

  How far either side to look, `1e-5` by default. Large enough that the
  one-sided derivatives are not dominated by rounding, small enough that
  a second feature does not fall inside the window.

## Value

A single logical, `TRUE` when the slope jumps.

## Details

The comparison is `abs(grad_y(at + eps) - grad_y(at - eps))` against
`1e-6 * max(1, abs(at))`, so the threshold is absolute near the origin
and relative away from it. At a genuine kink the difference is of the
order of the jump itself and does not shrink with `eps`: for a Laplace
at zero it is twice the rate. At a smooth point it is \\O(\varepsilon)\\
and falls well below the threshold.

A parent that cannot be differentiated at either point returns `NA` from
[`distributions7::distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.html)
or raises, and both are caught and read as `FALSE`: an unconfirmed
candidate is dropped rather than assumed.

## See also

[`distrib_kinks()`](https://statmodels7.github.io/penalties7/reference/distrib_kinks.md),
[`probe_theta()`](https://statmodels7.github.io/penalties7/reference/probe_theta.md)
