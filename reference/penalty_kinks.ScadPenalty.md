# Smoothness and Kind of a SCAD Penalty

[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
returns the five points where SCAD changes branch, and
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
returns `FALSE`, SCAD being bounded and so not a density at any
constant.

## Arguments

- pen:

  A
  [`ScadPenalty()`](https://statmodels7.github.io/penalties7/reference/ScadPenalty.md)
  object.

- theta:

  A list holding `lambda` and `a`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
a numeric vector of five points, in the order `0`, \\-\lambda\\,
\\\lambda\\, \\-a\lambda\\, \\a\lambda\\.
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
the single logical `FALSE`.

## Details

The five points are `0`, \\\pm\lambda\\ and \\\pm a\lambda\\, and they
are not alike. At the origin the first derivative jumps from
\\-\lambda\\ to \\\lambda\\, which is the kink that produces exact
zeros. At the other four the first derivative is continuous, both
branches agreeing, and the second jumps between \\0\\ and \\-1/(a-1)\\.
All five are returned because a numerical Hessian straddling any of them
is as wrong as a numerical gradient straddling the origin.

Four of the five move with the hyperparameters, so `theta` is not
optional.

## See also

[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
and
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
for the generics,
[`penalty_kinks.McpPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.McpPenalty.md)
for the sibling family's three points.

## Examples

``` r
penalty_kinks(scad_penalty(), list(lambda = 1, a = 3.7))
#> [1]  0.0 -1.0  1.0 -3.7  3.7
penalty_kinks(scad_penalty(), list(lambda = 2, a = 3.7))
#> [1]  0.0 -2.0  2.0 -7.4  7.4
is_proper(scad_penalty())
#> [1] FALSE
```
