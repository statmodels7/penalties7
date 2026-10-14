# Test for a Jump in the Slope of the Log-Density

Compares the derivative of the log-density in its argument on either
side of a candidate point at two widths, and returns `TRUE` when the
difference has the behavior of a jump.
[`distrib_kinks()`](https://statmodels7.github.io/penalties7/reference/distrib_kinks.md)
uses it to confirm a candidate kink inferred from `params_smooth`.

## Usage

``` r
has_jump(d, theta, at, eps = 1e-05)
```

## Arguments

- d:

  A distributions7 object.

- theta:

  A named list of parameter values at which the log-density is
  evaluated, as
  [`probe_theta()`](https://statmodels7.github.io/penalties7/reference/probe_theta.md)
  returns.

- at:

  The candidate point, a single number.

- eps:

  The relative narrow width, `1e-5` by default; the wide width is ten
  times as large, and both are multiplied by \\\max(1, \lvert
  a\rvert)\\. It should be large enough that the one-sided derivatives
  are not dominated by rounding, and small enough that a second feature
  of the density does not fall inside the window.

## Value

A single logical, `TRUE` when the slope jumps.

## Details

With \\s = \max(1, \lvert a\rvert)\\, where \\a\\ is `at`, the
difference \\g(h) = \lvert\ell^{(y)}(a + h) - \ell^{(y)}(a - h)\rvert\\
is computed at the narrow width \\h = \varepsilon s\\ and at the wide
width \\h = 10\varepsilon s\\. At a kink the difference is of the order
of the jump at both widths (twice the rate for a Laplace at its
location) and does not grow with the width. At a smooth point it is
about \\2h\lvert\ell^{(yy)}\rvert\\, which can exceed any fixed
threshold when the curvature is large, and it grows tenfold from the
narrow width to the wide one. A kink is therefore declared when
\\g(\varepsilon s)\\ exceeds \\10^{-6}s\\, a threshold that is absolute
near the origin and relative away from it, and is more than half of
\\g(10\varepsilon s)\\. The widths are scaled up from \\\varepsilon\\
and not down, so that a slope that the parent obtains by its own
numerical differentiation still resolves the jump at the narrow width.

When
[`distributions7::distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.html)
returns a value that is not finite or signals an error at any of the
four points, the result is `FALSE`, so the candidate is dropped.

## See also

[`distrib_kinks()`](https://statmodels7.github.io/penalties7/reference/distrib_kinks.md),
[`probe_theta()`](https://statmodels7.github.io/penalties7/reference/probe_theta.md)
