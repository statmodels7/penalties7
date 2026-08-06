# Construct the SCAD and MCP Penalties

SCAD (Fan and Li, 2001) is defined by its derivative on \\t \ge 0\\,
\$\$\rho'(t) = \lambda \\\\ (t \le \lambda), \qquad \rho'(t) =
\dfrac{a\lambda - t}{a - 1} \\\\ (\lambda \< t \le a\lambda), \qquad
\rho'(t) = 0 \\\\ (t \> a\lambda),\$\$ extended evenly, with \\a \> 2\\;
MCP (Zhang, 2010) by \\\rho'(t) = (\lambda - t/\gamma)\_+\\ with
\\\gamma \> 1\\. Every quantity – the value, the derivatives in
\\\beta\\ and in the hyperparameters, and the mixed block – is a closed
piecewise form. Both report kinks at zero, where the second derivative
additionally jumps at the region boundaries;
[`check_penalty`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
keeps its grids away from all of them by asking the object.

## Usage

``` r
scad_penalty(
  map = NULL,
  n_coef = 1L,
  link_lambda = linkfunctions7::log_link(),
  link_a = linkfunctions7::bounded_link(lwr = 2)
)

mcp_penalty(
  map = NULL,
  n_coef = 1L,
  link_lambda = linkfunctions7::log_link(),
  link_gamma = linkfunctions7::bounded_link(lwr = 1)
)
```

## Arguments

- map:

  The matrix \\D\\, or `NULL` (default) for the identity.

- n_coef:

  The number of coefficients when `map` is `NULL`.

- link_lambda:

  The link carrying \\\lambda\\.

- link_a, link_gamma:

  The link carrying the shoulder parameter, lower bounded at 2 (SCAD) or
  1 (MCP).

## Value

An object of class `ScadPenalty` or `McpPenalty`.

## References

Fan, J. and Li, R. (2001). Variable selection via nonconcave penalized
likelihood and its oracle properties. *JASA* 96, 1348-1360.

Zhang, C.-H. (2010). Nearly unbiased variable selection under minimax
concave penalty. *Annals of Statistics* 38, 894-942.

## Examples

``` r
pen <- scad_penalty(n_coef = 3)
penalty_value(pen, c(0.5, 2, 5), list(lambda = 1, a = 3.7))
#> [1] 4.664815
is_proper(pen)
#> [1] FALSE
```
