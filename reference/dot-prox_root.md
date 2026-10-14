# The Proximal Operator of a Separable Penalty, by Its Stationary Points

Returns, for every coordinate, the stationary point of smallest
objective \\(\beta - v_j)^2/(2t) - \log f(\beta)\\, where \\f\\ is the
parent's density, among the stationary points found on the segment
between \\v_j\\ and the mode of the parent. When the parent is unimodal,
this point is the global minimizer. The function serves the separable
penalties that have no closed-form operator, such as the heavy-tailed
Student t prior.

## Usage

``` r
.prox_root(pen, v, step, theta, n_grid = 64L)

.prox_mode(ly)

.prox_root_widened(pen, h, vi, m = 0)

.prox_uniroot(h, lo, hi, flo, fhi)
```

## Arguments

- pen:

  A
  [`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
  object with no map and a block width of one.

- v:

  A numeric vector, the points pulled towards.

- step:

  The step length \\t\\, a single positive finite number.

- theta:

  The aligned hyperparameter list.

- n_grid:

  The number of grid points on each segment, ends included.

- ly:

  The function \\\ell^{(y)}\\ of the parent, vectorized.

- h:

  The function \\h\\ of one coordinate.

- vi:

  That coordinate's point.

- m:

  The parent's mode, as `.prox_mode()` returns it.

- lo, hi:

  The ends of a bracket, in either order.

- flo, fhi:

  The values of `h` there, of opposite signs.

## Value

A numeric vector of the same length as `v`.

`.prox_mode()` returns a single number, the mode of the parent or zero.

## Details

The stationary points are the roots of \\h(\beta) = (\beta - v)/t -
\ell^{(y)}(\beta)\\, whose second term is
[`distributions7::distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.html).
For a density that is not log-concave \\h\\ can have more than one root
(up to three for the Student t, two of them local minima of the
objective), and the first root that a bracket encloses need not be the
global minimizer.

For a unimodal parent with mode \\m\\, the penalty is non-increasing up
to \\m\\ and non-decreasing after it, so the global minimizer lies on
the segment between \\v\\ and \\m\\. Moving a point beyond \\v\\ to
\\v\\, or a point beyond \\m\\ to \\m\\, lowers the quadratic term and
does not raise the penalty. The mode is the root of \\\ell^{(y)}\\,
found once per call because the parent and \\\theta\\ are the same for
every coordinate: when \\\ell^{(y)}(0)\\ is exactly zero the mode is
taken as zero, and otherwise the root is bracketed by doubling a step
from zero in the direction in which \\\ell^{(y)}\\ points, up to
\\2^{60}\\, and refined by
[`uniroot()`](https://rdrr.io/r/stats/uniroot.html) to full relative
precision. When no change of sign is found, or \\\ell^{(y)}\\ is not
finite, the mode is taken as zero and the segment is the one between
zero and \\v\\. The search moves away from zero only while
\\\ell^{(y)}\\ keeps its sign, so the root it finds is a local mode,
which is the mode when the parent is unimodal.

The segment is scanned on `n_grid` equally spaced points, both ends
included, with \\h\\ evaluated for every coordinate in one call. Every
interval over which \\h\\ changes sign is refined by
[`uniroot()`](https://rdrr.io/r/stats/uniroot.html), a grid point where
\\h\\ is exactly zero is kept as it is, and the candidate of smallest
objective is returned, with the log-density computed by
[`distributions7::distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.html).

When the segment shows no change of sign, the bracket starts one unit
beyond \\\min(v, m)\\ and \\\max(v, m)\\, each end that does not yet
have the required sign is moved out by the current width of the bracket
up to sixty times, and the single root found in the final bracket is
returned. For a parent with more than one mode, the result in either
case is a stationary point that need not be the global minimizer,
because the scan covers only the segment between \\v\\ and one local
mode.

The tolerance passed to
[`uniroot()`](https://rdrr.io/r/stats/uniroot.html) is `eps^0.75` times
the end of the bracket nearer to zero, or the smallest positive double
when the bracket has zero as an end or contains it.
[`uniroot()`](https://rdrr.io/r/stats/uniroot.html) adds `2 eps` times
the current iterate to that tolerance, so the root is found to a
relative accuracy whatever its magnitude.

## Errors

A message naming the penalty where neither the segment nor the widened
bracket shows a change of sign.

## See also

[`penalty_prox.DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.DistribPenalty.md)
