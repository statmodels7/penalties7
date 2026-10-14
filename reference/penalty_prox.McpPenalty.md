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

  The step length \\t\\, a single positive finite number with \\t c_j \<
  \gamma\\ in every coordinate (\\t \< \gamma\\ at unit curvature).

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

MCP has two pieces where SCAD has three, because its taper starts at the
origin. The first piece is computed as \\s\\(u -
t\lambda)\_{+}\\\gamma/(\gamma - t)\\, in which \\\gamma - t\\ is a
difference of exact operands, so that the value stays accurate as \\t\\
approaches \\\gamma\\.

The step enters through the division by \\1 - t/\gamma\\: the curvature
is the constant \\-1/\gamma\\ inside the taper, so the subproblem is
convex only while \\t \< \gamma\\. At or beyond that bound the operator
is set-valued, and the call signals an error that prints \\\gamma\\ and
the largest scaled step, defined below.

The formula is written for unit curvature. With the curvature \\c_j\\ of
[`mcp_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
(argument `curv`), the operator is the formula above with \\t c_j\\ in
place of \\t\\ and \\\lambda/c_j\\ in place of \\\lambda\\. Under a
diagonal map the step is \\t d_j^2\\, so the condition becomes \\t c_j
d_j^2 \< \gamma\\ in every coordinate, which is stricter than \\t \<
\gamma\\ where \\c_j d_j^2 \> 1\\ and looser where \\c_j d_j^2 \< 1\\.

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

# Inside the taper the soft-thresholded value is divided by 1 - t/gamma;
# beyond gamma * lambda the point is returned unchanged.
penalty_prox(mcp_penalty(n_coef = 3), c(1.5, 2.5, 6), 1, th)
#> [1] 0.75 2.25 6.00

# At the same lambda, MCP with gamma = 3 shrinks this point less than SCAD
# with a = 3.7.
penalty_prox(mcp_penalty(n_coef = 1), 2.5, 1, th)
#> [1] 2.25
penalty_prox(scad_penalty(n_coef = 1), 2.5, 1, list(lambda = 1, a = 3.7))
#> [1] 1.794118

# At the convexity bound the operator is set-valued, and the call signals
# an error.
try(penalty_prox(mcp_penalty(n_coef = 1), 2, 3, th))
#> Error : the MCP proximal operator needs step < gamma (3 here), with the step
#>   multiplied by the curvature and, under a diagonal map, by the squared
#>   diagonal entry; the largest scaled step is 3. Beyond that bound the
#>   subproblem is not convex and the operator is set-valued.
```
