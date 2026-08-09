# Construct a Structured Quadratic Penalty

The exact negative log-density of the Gaussian prior \\\beta \sim N(0,
\Omega(\theta)^{-})\\, with the precision \\\Omega\\ a parameters7
matrix parameter: \$\$\rho = \tfrac{1}{2}\\\beta'\Omega(\theta)\beta -
\tfrac{1}{2}\log\mathrm{pdet}\\\Omega(\theta) + \tfrac{r}{2}\log
2\pi.\$\$ This is the correlated prior
[`quadratic_penalty`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
cannot express: there the matrix is a constant and one scale multiplies
it, here the hyperparameters are the structure's free values and reach
every entry.

## Usage

``` r
structured_penalty(structure)
```

## Arguments

- structure:

  A parameters7 `matrix_parameter`, read as the PRECISION of the prior.
  A rank-deficient structure (an improper prior) is admitted;
  `is_proper` then answers `FALSE` and the constant uses the rank and
  the log pseudo-determinant.

## Value

An object of class `StructuredPenalty`.

## Details

The hyperparameters ARE the structure's free vector, which is
unconstrained by construction, so every link is the identity – the same
flattening convention the multivariate families of distributions7
follow, and for the same reason: the constraint lives in the structure,
where a scalar link cannot express it. Every derivative comes from the
structure's own contract: the theta gradient is
\\\tfrac{1}{2}\beta'A_k\beta - \tfrac{1}{2}\partial_k\log\mathrm{pdet}\\
with \\A_k\\ the structure's `param_d1`, the theta Hessian adds
`param_d2`, and the mixed block is \\A_k\beta\\.

There is no `map` argument, deliberately: a linear image of a structured
precision is a different precision, and composing it into the structure
– where its log-determinant stays exact – is the structure's business,
not this constructor's.

## See also

[`quadratic_penalty`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md),
[`additive_penalty`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md),
[`distrib_penalty`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md)

## Examples

``` r
# an AR(1) prior on four coefficients: three hyperparameters reach every
# entry of the precision
pen <- structured_penalty(parameters7::ar1(4, role = "precision"))
theta <- list(log_scale = 0.2, z_rho = 0.5)
penalty_value(pen, c(0.3, -0.1, 0.4, 0.2), theta)
#> [1] 3.858268
penalty_rank(pen)
#> [1] 4
```
