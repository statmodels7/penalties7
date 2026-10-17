# Is a Penalty a Proper Prior?

Returns `TRUE` when \\\exp(-\rho)\\ integrates to one over the
coefficients, so the value returned by
[`penalty_value()`](https://statmodels7.github.io/penalties7/reference/penalty_value.md)
is exactly the negative log-density of a prior in \\\beta\\ and may be
added to a negative log-likelihood as it stands. Returns `FALSE` when it
does not. In that case the value keeps the constant that its branch
defines, if any: a rank-deficient quadratic, additive or structured
penalty keeps its pseudo-determinant term, a separable penalty keeps the
parent's own constant without a Jacobian term, and SCAD and MCP carry no
constant.

## Usage

``` r
is_proper(pen, ...)
```

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object of any branch.

- ...:

  Passed to methods. The shipped methods do not read it.

## Value

A single logical.

## Details

Three situations make a penalty improper, and they arise from different
causes.

A rank-deficient quadratic puts no cost at all on its null directions,
so \\\exp(-\rho)\\ is flat along them and its integral diverges. This is
the ordinary case for a smoothing penalty: second differences over \\q\\
coefficients leave the level and the slope unpenalized and have rank
\\q - 2\\. The penalty is still usable, because the likelihood supplies
the missing curvature, and its normalizing constant is the log
pseudo-determinant over the range, which
[`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
returns. The same holds for a rank-deficient additive or structured
penalty.

A separable penalty is a density in \\t = D\beta\\. It is proper without
a map, and under a square nonsingular map, where the value carries the
Jacobian term \\-\log\lvert\det D\rvert\\ and is a density in \\\beta\\.
Under a map with fewer rows than columns, or a singular square map (a
diagonal map with a zero entry, for instance), \\\exp(-\rho)\\ is flat
along the null space of the map, and the predicate is `FALSE`.

SCAD and MCP are improper for a third reason. Both are defined by their
derivative, which is exactly zero beyond \\a\lambda\\ and
\\\gamma\lambda\\, so \\\rho\\ is constant in every direction far from
the origin and \\\exp(-\rho)\\ is not integrable for any constant. They
are not the negative log-density of any prior, and the quantities of a
marginal criterion are not available for them.

The predicate is a property of the penalty as constructed, so it does
not depend on `theta` and does not take it.

## See also

[`penalty_value()`](https://statmodels7.github.io/penalties7/reference/penalty_value.md)
for the value whose constant this predicate describes,
[`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
for the pseudo-determinant of a deficient quadratic,
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
and
[`has_prox()`](https://statmodels7.github.io/penalties7/reference/has_prox.md)
for the other two predicates on which a fitting routine selects a
method.

## Examples

``` r
# A full-rank quadratic is a proper Gaussian prior.
is_proper(quadratic_penalty(diag(2)))
#> [1] TRUE

# Second differences leave the level and the slope unpenalized.
is_proper(quadratic_penalty(crossprod(diff(diag(4), differences = 2))))
#> [1] FALSE
penalty_rank(quadratic_penalty(crossprod(diff(diag(4), differences = 2))))
#> [1] 2

# The lasso is a proper Laplace prior; SCAD and MCP are not densities.
c(lasso = is_proper(lasso_penalty()),
  scad  = is_proper(scad_penalty()),
  mcp   = is_proper(mcp_penalty()))
#> lasso  scad   mcp 
#>  TRUE FALSE FALSE 

# A lasso on second differences leaves a two-dimensional null space.
is_proper(lasso_penalty(map = diff(diag(4), differences = 2)))
#> [1] FALSE
```
