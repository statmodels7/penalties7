# Proximal Operator of SCAD

The closed piecewise operator of Fan and Li (2001): a soft threshold
near zero, a rescaled threshold over the tapering region, and the
identity beyond \\a\lambda\\, so a large coefficient passes through
untouched.

## Arguments

- pen:

  A
  [`ScadPenalty()`](https://statmodels7.github.io/penalties7/reference/ScadPenalty.md)
  object.

- v:

  A numeric vector of length `pen@n_coef`.

- step:

  The step length \\t\\, a single positive number below \\a-1\\.

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

The middle piece divides by \\1 - t/(a-1)\\, which is where the step
matters: SCAD's curvature is \\-1/(a-1)\\ there, so the subproblem is
convex only while \\t \< a-1\\. At or beyond that the operator is
set-valued and the call is rejected with the bound printed. Under a
diagonal map the condition tightens to \\t \< (a-1)/d_j^2\\.

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

# One point in each region: thresholded, rescaled, untouched.
penalty_prox(scad_penalty(n_coef = 3), c(1.5, 3, 6), 1, th)
#> [1] 0.500000 2.588235 6.000000

# Past a * lambda the operator is the identity, so a large coefficient is
# not shrunk at all. The lasso shrinks the same point by t * lambda.
penalty_prox(scad_penalty(n_coef = 1), 6, 1, th)
#> [1] 6
penalty_prox(lasso_penalty(n_coef = 1), 6, 1, list(lambda = 1))
#> [1] 5

# Beyond the convexity bound the operator is set-valued and rejects.
try(penalty_prox(scad_penalty(n_coef = 1), 2, 2.7, th))
#> Error : the SCAD proximal operator needs step < a - 1 (2.7 here): beyond that the
#>   subproblem is not convex and the operator is set-valued.
```
