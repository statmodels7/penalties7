# Value of a SCAD Penalty

Returns the closed antiderivative of SCAD's defining \\\rho'\\, summed
over the coordinates of \\D\beta\\ and anchored at \\\rho(0) = 0\\. No
constant is added, SCAD being improper.

## Arguments

- pen:

  A
  [`ScadPenalty()`](https://statmodels7.github.io/penalties7/reference/ScadPenalty.md)
  object.

- beta:

  A numeric vector of length `pen@n_coef`, already coerced by the
  generic.

- theta:

  A list holding `lambda` and `a`, already aligned and bound-checked by
  the generic.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A single number, the sum over coordinates. Bounded above by
`pen@n_coef * (a + 1) * lambda^2 / 2`.

## Details

With \\u = \lvert t\rvert\\ on each coordinate,

\$\$\rho(u) = \lambda u \\\\ (u \le \lambda), \qquad \rho(u) =
\frac{2a\lambda u - u^2 - \lambda^2}{2(a-1)} \\\\ (\lambda \< u \le
a\lambda), \qquad \rho(u) = \frac{(a+1)\lambda^2}{2} \\\\ (u \>
a\lambda),\$\$

the three pieces meeting continuously at \\\lambda\\ and \\a\lambda\\.
The value saturates, leaving a large coefficient unshrunk and making the
penalty improper.

## See also

[`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
for the construction,
[`penalty_gradient.ScadPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.ScadPenalty.md)
for the derivative it integrates,
[`penalty_value.McpPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.McpPenalty.md)
for the sibling family.

## Examples

``` r
pen <- scad_penalty(n_coef = 1)
th <- list(lambda = 1, a = 3.7)

# One coordinate in each region: linear, tapering, flat.
sapply(c(0.5, 2, 100), function(t) penalty_value(pen, t, th))
#> [1] 0.500000 1.814815 2.350000

# The flat value is (a + 1) lambda^2 / 2.
(3.7 + 1) / 2
#> [1] 2.35
```
