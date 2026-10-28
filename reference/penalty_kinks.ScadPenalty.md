# Smoothness and Kind of a SCAD Penalty

[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
returns the points where SCAD changes branch, and
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
returns `FALSE`, since SCAD is bounded and so not a density at any
constant.

## Arguments

- penalty:

  A
  [`ScadPenalty()`](https://statmodels7.github.io/penalties7/reference/ScadPenalty.md)
  object.

- theta:

  A list holding `lambda` and `a`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
returns a numeric vector of distinct points: `0`, then \\-\lambda/c_j\\
and \\\lambda/c_j\\ over the distinct rates, then \\-a\lambda/c_j\\ and
\\a\lambda/c_j\\ over the same rates. At unit curvature these are the
five points `0`, \\-\lambda\\, \\\lambda\\, \\-a\lambda\\, \\a\lambda\\,
in this order.
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
returns the single logical `FALSE`.

## Details

At unit curvature the points are `0`, \\\pm\lambda\\ and \\\pm
a\lambda\\, and they are of two kinds. At the origin the first
derivative jumps from \\-\lambda\\ to \\\lambda\\, which is the kink
that produces exact zeros. At the other four the first derivative is
continuous and the second jumps between \\0\\ and \\-1/(a-1)\\
(\\-c_j/(a-1)\\ with curvature \\c_j\\). All of them are returned
because a numerical Hessian across any of them is as wrong as a
numerical gradient across the origin.

With `curvature`, the points \\\pm\lambda/c_j\\ and \\\pm a\lambda/c_j\\
are returned for each distinct value of \\\lambda/c_j\\, so the vector
has more than five points when the curvatures differ. Every point except
the origin moves with the hyperparameters, so `theta` is required.

## See also

[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
and
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
for the generics,
[`penalty_kinks.McpPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.McpPenalty.md)
for the sibling family.

## Examples

``` r
penalty_kinks(scad_penalty(), list(lambda = 1, a = 3.7))
#> [1]  0.0 -1.0  1.0 -3.7  3.7
penalty_kinks(scad_penalty(), list(lambda = 2, a = 3.7))
#> [1]  0.0 -2.0  2.0 -7.4  7.4

# Two curvatures give two rates and nine points.
penalty_kinks(scad_penalty(n_coef = 2, curvature = c(1, 2)),
              list(lambda = 1, a = 3.7))
#> [1]  0.00 -1.00 -0.50  1.00  0.50 -3.70 -1.85  3.70  1.85
is_proper(scad_penalty())
#> [1] FALSE
```
