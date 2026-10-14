# Marginal Quantities of a Structured Penalty

The four quantities that a REML or marginal likelihood criterion reads.
[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
returns the precision \\\Omega(\theta)\\,
[`penalty_rank()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
the structure's rank,
[`penalty_null_basis()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
its null basis, and
[`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
the log pseudo-determinant of the precision with its first two
derivatives in the free values.

## Arguments

- pen:

  A
  [`StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/StructuredPenalty.md)
  object.

- theta:

  A named list of the structure's free values.
  [`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
  and
  [`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
  only.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
the symmetric precision matrix, of side `pen@n_coef`, carrying the
structure's own dimnames.
[`penalty_rank()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
a single integer, the structure's rank.
[`penalty_null_basis()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
a matrix with `pen@n_coef` rows and one column per null direction, with
no columns for a full-rank structure.
[`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
a list of `value` (a number), `grad` (a list of one number per free
value) and `hess` (a list of one number per unordered pair, keyed
diagonals first).

## Details

For a plain quadratic penalty \\\log^{+}\lvert S\rvert\\ is linear in
\\\log\lambda\\. Here the log-determinant depends on every free value,
and its value and derivatives are those of the structure's
`param_logdet()`, `param_dlogdet()` and `param_d2logdet()`.

The rank and the null basis are properties of the structure and do not
depend on \\\theta\\: a `matrix_parameter` records them at construction,
so a family whose rank is deficient is deficient in the same directions
at every free vector.
[`penalty_null_basis()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
reshapes a vector-valued null basis into a one-column matrix, so its
return is always a matrix.

## See also

[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
and its three siblings for the generics,
[`struct_omega()`](https://statmodels7.github.io/penalties7/reference/struct_omega.md)
for where the derivatives come from,
[`penalty_matrix.QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.QuadraticPenalty.md)
for the branch whose log pseudo-determinant is linear in
\\\log\lambda\\.

## Examples

``` r
pen <- structured_penalty(parameters7::ar1(4))
th <- list(log_scale = 0.2, z_rho = 0.5)

penalty_rank(pen)
#> [1] 4
dim(penalty_null_basis(pen))
#> [1] 4 0

# The matrix is the precision itself, so it is the Hessian of the value.
all.equal(penalty_matrix(pen, th), penalty_hessian(pen, c(0, 0, 0, 0), th))
#> [1] TRUE

# And its log pseudo-determinant is the ordinary log-determinant here,
# the structure being full rank.
lp <- penalty_logpdet(pen, th)
lp$value - as.numeric(
  determinant(penalty_matrix(pen, th), logarithm = TRUE)$modulus)
#> [1] -3.330669e-16
names(lp$hess)
#> [1] "log_scale_log_scale" "z_rho_z_rho"         "log_scale_z_rho"    
```
