# Proximal Operator of MCP

The closed piecewise operator of Zhang (2010): a rescaled soft threshold
below \\\gamma\lambda\\ and the identity beyond it, so a large
coefficient passes through untouched.

## Arguments

- pen:

  An
  [`McpPenalty()`](https://statmodels7.github.io/penalties7/reference/ScadPenalty.md)
  object.

- v:

  A numeric vector of length `pen@n_coef`.

- step:

  The step length \\t\\, a single positive number below \\\gamma\\.

- theta:

  A named list holding `lambda` and `gamma`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of the same length as `v`.

## Details

With \\u = \lvert v\rvert\\ and \\s\\ its sign,

\$\$\mathrm{prox}(v) = \begin{cases} \dfrac{s\\(u - t\lambda)\_{+}}{1 -
t/\gamma} & u \le \gamma\lambda \\ v & u \> \gamma\lambda.
\end{cases}\$\$

MCP has two pieces where SCAD has three, its taper starting at the
origin. The division by \\1 - t/\gamma\\ is where the step matters: the
curvature is a constant \\-1/\gamma\\ inside the taper, so the
subproblem is convex only while \\t \< \gamma\\. At or beyond that the
operator is set-valued and the call is rejected with the bound printed.
Under a diagonal map the condition tightens to \\t \< \gamma/d_j^2\\.

## References

Zhang, C.-H. (2010). Nearly unbiased variable selection under minimax
concave penalty. *Annals of Statistics* **38**, 894-942.

## See also

[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
for the generic,
[`mcp_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
for the penalty,
[`penalty_prox.ScadPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.ScadPenalty.md)
for the sibling family,
[`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
for the same operator as a table.

## Examples

``` r
th <- list(lambda = 1, gamma = 3)

# Inside the taper the threshold is rescaled upward by 1/(1 - t/gamma);
# beyond gamma * lambda the point is returned unchanged.
penalty_prox(mcp_penalty(n_coef = 3), c(1.5, 2.5, 6), 1, th)
#> [1] 0.75 2.25 6.00

# The rescaling lets MCP reach the unpenalized answer sooner than SCAD at
# the same lambda.
penalty_prox(mcp_penalty(n_coef = 1), 2.5, 1, th)
#> [1] 2.25
penalty_prox(scad_penalty(n_coef = 1), 2.5, 1, list(lambda = 1, a = 3.7))
#> [1] 1.794118

# Beyond the convexity bound the operator is set-valued and rejects.
try(penalty_prox(mcp_penalty(n_coef = 1), 2, 3, th))
#> Error : the MCP proximal operator needs step < gamma (3 here): beyond that the
#>   subproblem is not convex and the operator is set-valued.
```
