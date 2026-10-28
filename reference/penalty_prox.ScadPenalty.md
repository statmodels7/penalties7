# Proximal Operator of SCAD

Computes the closed piecewise operator of Fan and Li (2001): a soft
threshold near zero, a rescaled threshold over the tapering region, and
the identity beyond \\a\lambda\\, so a large coefficient passes through
untouched.

## Arguments

- penalty:

  A
  [`ScadPenalty()`](https://statmodels7.github.io/penalties7/reference/ScadPenalty.md)
  object.

- v:

  A numeric vector of length `penalty@n_coef`.

- step:

  The step length \\t\\, a single positive finite number with \\t c_j \<
  a-1\\ in every coordinate (\\t \< a - 1\\ at unit curvature).

- theta:

  A named list holding `lambda` and `a`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of the same length as `v`.

## Details

With \\u = \lvert v\rvert\\ and \\s\\ its sign,

\$\$\mathrm{prox}(v) = \begin{cases} s\\(u - t\lambda)\_{+} & u \le
(1+t)\lambda \\ \dfrac{s\\(u - t a \lambda/(a-1))}{1 - t/(a-1)} &
(1+t)\lambda \< u \le a\lambda \\ v & u \> a\lambda. \end{cases}\$\$

The middle piece is computed in the equivalent form \\s\\(u -
t(a\lambda - u)/(a - 1 - t))\\, in which \\a\lambda - u\\ and \\a - 1 -
t\\ are differences of exact operands, so that the value stays accurate
as \\t\\ approaches \\a - 1\\.

The step enters through the division by \\1 - t/(a-1)\\: the curvature
of SCAD is \\-1/(a-1)\\ in the middle piece, so the subproblem is convex
only while \\t \< a-1\\. At or beyond that bound the operator is
set-valued, and the call signals an error that prints \\a - 1\\ and the
largest scaled step, defined below.

The formula is written for unit curvature. With the curvature \\c_j\\ of
[`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
(argument `curvature`), the penalty in coordinate \\j\\ is
\\c_j\rho(\cdot;\\ \lambda/c_j)\\, and the operator is the formula above
with \\t c_j\\ in place of \\t\\ and \\\lambda/c_j\\ in place of
\\\lambda\\. Under a diagonal map the step is \\t d_j^2\\, so the
condition becomes \\t c_j d_j^2 \< a - 1\\ in every coordinate, which is
stricter than \\t \< a - 1\\ where \\c_j d_j^2 \> 1\\ and looser where
\\c_j d_j^2 \< 1\\.

## References

Fan, J. and Li, R. (2001). Variable selection via nonconcave penalized
likelihood and its oracle properties. *Journal of the American
Statistical Association* **96**, 1348-1360.

## See also

[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
for the generic,
[`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
for the penalty,
[`penalty_prox.McpPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.McpPenalty.md)
for the sibling family,
[`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
for the same operator as a table.

## Examples

``` r
th <- list(lambda = 1, a = 3.7)

# The three points lie one in each region: the first is thresholded, the
# second rescaled and the third untouched.
penalty_prox(scad_penalty(n_coef = 3), c(1.5, 3, 6), 1, th)
#> [1] 0.500000 2.588235 6.000000

# Past a * lambda the operator is the identity, so a large coefficient is
# not shrunk at all. The lasso shrinks the same point by t * lambda.
penalty_prox(scad_penalty(n_coef = 1), 6, 1, th)
#> [1] 6
penalty_prox(lasso_penalty(n_coef = 1), 6, 1, list(lambda = 1))
#> [1] 5

# At the convexity bound the operator is set-valued, and the call signals
# an error.
try(penalty_prox(scad_penalty(n_coef = 1), 2, 2.7, th))
#> Error : the SCAD proximal operator needs step < a - 1 (2.7 here), where the
#>   step is multiplied by the curvature and, under a diagonal map, by
#>   the squared diagonal entry; the largest scaled step is 2.7. Beyond
#>   that bound the subproblem is not convex and the operator is
#>   set-valued.
```
