# Marginal Derivatives of a Structured Penalty

The methods of the structured branch for the four generics that a
marginal criterion requires beyond the second order. \\S =
\Omega(\theta)\\ is the matrix parameter itself, so
[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
and
[`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md)
are the derivatives `param_d1` and `param_d2` of the structure,
[`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md)
is the second of those applied to the coefficients, and
[`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md)
is `TRUE`.

## Arguments

- pen:

  A
  [`StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/StructuredPenalty.md)
  object.

- beta:

  A numeric vector of length `pen@n_coef`. Read by
  [`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md);
  the other three methods ignore it.

- theta:

  A named list of the free values of the structure.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
a list of one matrix per free value, named by `pen@params` and with
dimnames stripped.
[`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md)
a list of one matrix per unordered pair, keyed diagonals first.
[`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md)
a list of one numeric vector of length `pen@n_coef` per unordered pair,
keyed the same way.
[`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md)
the single logical `TRUE`.

## Details

Nothing is differentiated here. parameters7 supplies both derivative
arrays in closed form for the structures that it ships, and through its
numerical fallback for a structure that does not implement them; this
branch unwraps and re-keys the arrays. The free values in `theta` are
passed to the structure as they are, in `pen@params` order.

Since the value is quadratic in the coefficients, the mixed third
derivative is
\\(\partial^2\Omega/\partial\theta_m\partial\theta_l)\beta\\, so
[`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md)
calls
[`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md)
and multiplies.

## The keying

`param_d2` is a **flat** list whose keys join the free names with a
separator of its own, so an entry is located through
[`parameters7::param_tuple_indices()`](https://statmodels7.github.io/parameters7/reference/param_tuple_indices.html),
which enumerates the tuples in the same order as the components. The
method does not parse the keys, because a free name may itself contain
the separator. If a pair has no component, the method signals an error.

## See also

[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
and its siblings for the generics,
[`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md)
for the branch,
[`parameters7::param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.html)
for the arrays that these methods read,
[`penalty_dhessian_beta.StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian_beta.penalty.md)
and
[`penalty_d2hessian_beta.StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.penalty.md)
for the derivatives of the Hessian in the coefficients.

## Examples

``` r
s <- structured_penalty(parameters7::log_cholesky(2))
th <- list(log_L1 = 0.1, log_L2 = -0.1, L2.1 = 0.3)
b <- c(1, -0.5)

# Three free values, so three first derivatives and six pairs.
names(penalty_dhessian(s, b, th))
#> [1] "log_L1" "log_L2" "L2.1"  
names(penalty_d2hessian(s, b, th))
#> [1] "log_L1_log_L1" "log_L2_log_L2" "L2.1_L2.1"     "log_L1_log_L2"
#> [5] "log_L1_L2.1"   "log_L2_L2.1"  

# Unlike the quadratic branch, the second order does not vanish.
penalty_d2hessian(s, b, th)$log_L1_log_L1
#>           [,1]      [,2]
#> [1,] 4.8856110 0.3315513
#> [2,] 0.3315513 0.0000000

# And the mixed block is that matrix applied to the coefficients.
all.equal(penalty_dcross(s, b, th)$log_L1_log_L1,
          drop(penalty_d2hessian(s, b, th)$log_L1_log_L1 %*% b))
#> [1] TRUE
```
