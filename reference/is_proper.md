# Is a Penalty a Proper Prior?

`TRUE` when \\\exp(-\rho)\\ integrates to one over the coefficients, so
the value returned by
[`penalty_value()`](https://statmodels7.github.io/penalties7/reference/penalty_value.md)
is exactly the negative log-density of a prior and may be added to a
negative log-likelihood as it stands. `FALSE` when it does not, and then
the value is the bare \\\rho\\ with no constant attached.

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

  Passed to methods. No shipped method reads it.

## Value

A single logical.

## Details

Two things make a penalty improper, and they are different failures.

A **rank-deficient quadratic** puts no cost at all on its null
directions, so \\\exp(-\rho)\\ is flat along them and its integral
diverges. This is the ordinary case for a smoothing penalty: second
differences over \\q\\ coefficients leave the level and the slope free
and have rank \\q - 2\\. The penalty is still usable, because the
likelihood supplies the missing curvature; what changes is that the
normalizing constant is the log **pseudo**-determinant over the range,
which
[`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
returns.

**SCAD and MCP** are improper for a different reason: both are defined
by their derivative, which is exactly zero beyond \\a\lambda\\ and
\\\gamma\lambda\\, so \\\rho\\ is flat in every direction far from the
origin and no constant makes \\\exp(-\rho)\\ integrable. There is no
density behind them at all, and no marginal criterion reaches them.

The predicate is a statement about the penalty as constructed, so it
does not depend on `theta` and takes none.

## See also

[`penalty_value()`](https://statmodels7.github.io/penalties7/reference/penalty_value.md)
for the constant this is about,
[`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
for the pseudo-determinant a deficient quadratic uses,
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
and
[`has_prox()`](https://statmodels7.github.io/penalties7/reference/has_prox.md)
for the other two predicates a consumer routes on.

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

# The lasso is a proper Laplace prior; SCAD and MCP are densities of
# nothing.
c(lasso = is_proper(lasso_penalty()),
  scad  = is_proper(scad_penalty()),
  mcp   = is_proper(mcp_penalty()))
#> lasso  scad   mcp 
#>  TRUE FALSE FALSE 
```
