# Construct a Separable Penalty From a Distribution

\\\rho(D\beta;\theta) = -\sum_i \log f(b_i;\theta)\\ for a
distributions7 density \\f\\ read at the successive blocks \\b_i\\ of
\\D\beta\\. A univariate parent gives blocks of one coordinate, which is
the separable penalty; a \\p\\-variate parent gives blocks of \\p\\,
which is a prior that lets the coordinates of one block depend on each
other while the blocks stay independent.

The hyperparameters ARE the distribution's free parameters; bounds and
links are read off the distribution object rather than restated, and
every derivative is the distribution's, reassembled: the gradient is
\\-D'\ell^{(y)}\\, the Hessian \\-D'\mathrm{diag}(\ell^{(yy)})D\\ (block
diagonal rather than diagonal when the parent is multivariate), the
theta blocks the summed score and Hessian, and the mixed block one
\\-D'\ell^{(y\theta_k)}\\ per hyperparameter – which is what
[`distrib_cross_y`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.html)
exists for.

## Usage

``` r
distrib_penalty(d, map = NULL, n_coef = NULL, kinks = NULL)
```

## Arguments

- d:

  A continuous distributions7 object; typically a `fixed()` wrapper
  holding the location at zero. A multivariate parent of dimension \\p\\
  is read blockwise, and the number of coefficients must then be a
  multiple of \\p\\.

- map:

  The matrix \\D\\, or `NULL` (default) for the identity.

- n_coef:

  The number of coefficients; required when `map` is `NULL`, ignored
  otherwise.

- kinks:

  The points where the parent's log-density is not differentiable in its
  argument. `NULL`, the default, derives them from the parent with
  [`distrib_kinks`](https://statmodels7.github.io/penalties7/reference/distrib_kinks.md);
  pass a numeric vector to say so directly, or `numeric(0)` to declare
  there are none. A multivariate parent has none: a kink is a point of a
  scalar argument.

## Value

An object of class `DistribPenalty`.

## Details

Lasso, the elastic net and the heavy-tailed prior are this construction
at a
[`fixed`](https://statmodels7.github.io/distributions7/reference/fixed.html)
Laplace, elastic net and Student t; see
[`ridge_penalty`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md).
The normalizing constant comes with the density and is kept, so the
value is exactly the negative log prior density and a free scale or a
free \\\nu\\ is estimable.

A multivariate parent is centered by the caller, typically through
[`fixed`](https://statmodels7.github.io/distributions7/reference/fixed.html)
at a zero mean, and its matrix parameter carries the dependence within a
block. It has no proximal operator: that operator acts one coordinate at
a time and the coordinates of a block do not separate.

## See also

[`quadratic_penalty`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md),
[`additive_penalty`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md),
[`structured_penalty`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md)

## Examples

``` r
d <- distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0)
pen <- distrib_penalty(d, n_coef = 3)
penalty_value(pen, c(1, 0, -1), list(sigma = 2))
#> [1] 5.086257

# a correlated prior over two coefficients per block, three blocks
mv <- distributions7::fixed(distributions7::mvgaussian_distrib(2),
                            mu1 = 0, mu2 = 0)
pen2 <- distrib_penalty(mv, n_coef = 6)
penalty_value(pen2, c(1, 0, -1, 0.5, 0.2, -0.3),
              list(sigma_log_L1 = 0, sigma_log_L2 = 0, sigma_L2.1 = 0.4))
#> [1] 7.090831
```
