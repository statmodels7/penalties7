# S7 Class for Smoothers of the Absolute Value

A smooth replacement \\s(u)\\ for \\\lvert u \rvert\\, carrying its
derivatives in \\u\\ up to order five as functions. Every non-smooth
primitive the toolkit uses is generated from the absolute value –

\$\$\operatorname{sign}(u) = \frac{d\lvert u\rvert}{du}, \qquad
\mathbb{1}(u \geq 0) = \frac{1 + \operatorname{sign}(u)}{2}, \qquad
(u)\_{+} = \frac{u + \lvert u\rvert}{2},\$\$

so one contract serves them all: the smooth sign is \\s'\\, the smooth
step \\(1 + s'(u))/2\\ and the smooth hinge \\(u + s(u))/2\\ follow by
composition, and a consumer that replaces \\\lvert u\rvert\\ by \\s(u)\\
has replaced every one of them consistently. A break-point term smoothed
this way becomes an ordinary nonlinear term whose block is the true
Jacobian, which is what makes a random or penalized development of its
break-points fittable; a kinked penalty smoothed this way becomes a
proper separable one, at the price of the exact zeros the kink produced.

## Usage

``` r
abs_smoother(
  smoother_name = character(0),
  width = NULL,
  width_name = "h",
  per_group = FALSE,
  s = list(),
  width_from_spacing = NULL,
  tau_correction = NULL,
  exact_radius = NULL
)
```

## Arguments

- smoother_name:

  A string naming the smoother.

- width:

  The transition width, or `NULL` to be resolved at build.

- width_name:

  What the width parameter is called (`"h"` or `"c"`).

- per_group:

  Whether the width is resolved per group.

- s:

  A list of six functions of `(u, width)`: \\s\\ and its derivatives in
  \\u\\ of orders one to five, in order.

- width_from_spacing:

  `NULL` (the identity) or a function carrying a spacing in covariate
  units onto the width parameter's own scale (the square, for the
  hyperbolic).

- tau_correction:

  `NULL`, or a function `(tau, width)` returning the corrected scale of
  a random break-point.

- exact_radius:

  `NULL`, or a function of the width returning the radius beyond which
  \\s(u) = \lvert u\rvert\\ exactly.

## Value

An object of class `abs_smoother`.

## Details

The derivatives are FUNCTIONS and not expressions: a piecewise smoother
(the quintic) has branches, which
[`deriv`](https://rdrr.io/r/stats/deriv.html) does not read, and with
the derivatives written in the constructor the branches are ordinary
code. Each function takes `(u, width)` and vectorizes in both, so a
per-group width is one value per observation.

`width` is the transition scale – `h` for a smoother whose parameter is
a length, the bent-cable reading of a transition of width \\h\\; `c` for
the hyperbolic, whose parameter is a squared length. `NULL`, the
default, asks the consumer to resolve it from the data at build (a
break-point term takes the median spacing of its covariate), through
[`smoother_width`](https://statmodels7.github.io/penalties7/reference/smoother_width.md).
`per_group` asks for one width per group where a grouping is available,
the validity window of a Laplace approximation being per-subject.

`tau_correction` is a property of the mollifier, not of any model: the
probit smoother satisfies an exact convolution identity, smoothing with
width \\h\\ being the same as convolving the break-point with \\N(0,
h^2)\\, so an apparent scale \\\tau\\ of a random break-point composes
as \\\tau^2\_{\mathrm{true}} = \tau^2 - h^2\\ and the smoother declares
the correction. Smoothers with no such identity declare `NULL` and a
consumer reports the apparent scale alone.

## References

Bacon, D. W. and Watts, D. G. (1971). Estimating the transition between
two intersecting straight lines. *Biometrika*, 58(3), 525–534.

Tishler, A. and Zang, I. (1981). A new maximum likelihood algorithm for
piecewise regression. *Journal of the American Statistical Association*,
76(376), 980–987.

Seo, M. H. and Linton, O. (2007). A smoothed least squares estimator for
threshold regression models. *Journal of Econometrics*, 141(2), 704–735.

## See also

[`smooth_probit`](https://statmodels7.github.io/penalties7/reference/smooth_probit.md),
[`smooth_hyperbolic`](https://statmodels7.github.io/penalties7/reference/smooth_hyperbolic.md),
[`smooth_quintic`](https://statmodels7.github.io/penalties7/reference/smooth_quintic.md),
[`check_abs_smoother`](https://statmodels7.github.io/penalties7/reference/check_abs_smoother.md)

## Examples

``` r
S7::S7_inherits(smooth_probit(), abs_smoother)
#> [1] TRUE
```
