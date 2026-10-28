# Value of a SCAD Penalty

Returns the closed antiderivative of SCAD's defining \\\rho'\\, summed
over the coordinates of \\D\beta\\ and anchored at \\\rho(0) = 0\\. It
adds no constant, since SCAD is improper.

## Arguments

- penalty:

  A
  [`ScadPenalty()`](https://statmodels7.github.io/penalties7/reference/ScadPenalty.md)
  object.

- beta:

  A numeric vector of length `penalty@n_coef`, already coerced by the
  generic.

- theta:

  A list holding `lambda` and `a`, already aligned and bound-checked by
  the generic.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A single number, the sum over coordinates. It is bounded above by
\\\sum_j (a+1)\lambda^2/(2c_j)\\, where the sum runs over the
coordinates of \\D\beta\\; under the identity map at unit curvature the
bound is `penalty@n_coef * (a + 1) * lambda^2 / 2`.

## Details

With unit curvature and \\u = \lvert t\rvert\\ on each coordinate,

\$\$\rho(u) = \lambda u \\\\ (u \le \lambda), \qquad \rho(u) =
\frac{2a\lambda u - u^2 - \lambda^2}{2(a-1)} \\\\ (\lambda \< u \le
a\lambda), \qquad \rho(u) = \frac{(a+1)\lambda^2}{2} \\\\ (u \>
a\lambda),\$\$

where the three pieces meet continuously at \\\lambda\\ and
\\a\lambda\\. With curvature \\c_j\\ (the `curvature` property),
coordinate \\j\\ contributes \\c_j\\\rho(u_j;\\ \lambda/c_j)\\, the same
formula with \\\lambda/c_j\\ in place of \\\lambda\\ and multiplied by
\\c_j\\. The value saturates, so a large coefficient receives no
shrinkage and the penalty is improper.

## See also

[`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
for the construction,
[`penalty_gradient.ScadPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.ScadPenalty.md)
for the derivative that it integrates,
[`penalty_value.McpPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.McpPenalty.md)
for the sibling family.

## Examples

``` r
pen <- scad_penalty(n_coef = 1)
th <- list(lambda = 1, a = 3.7)

# One coordinate is evaluated in each region: linear, tapering, flat.
sapply(c(0.5, 2, 100), function(t) penalty_value(pen, t, th))
#> [1] 0.500000 1.814815 2.350000

# The flat value is (a + 1) lambda^2 / 2.
(3.7 + 1) / 2
#> [1] 2.35
```
