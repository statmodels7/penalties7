# Construct a Separable Penalty From a Distribution

\\\rho(D\beta;\theta) = -\sum_j \log f((D\beta)\_j;\theta)\\ for a
univariate distributions7 density \\f\\. The hyperparameters ARE the
distribution's free parameters; bounds and links are read off the
distribution object rather than restated, and every derivative is the
distribution's, reassembled: the gradient is \\-D'\ell^{(y)}\\, the
Hessian \\-D'\mathrm{diag}(\ell^{(yy)})D\\, the theta blocks the summed
score and Hessian, and the mixed block one \\-D'\ell^{(y\theta_k)}\\ per
hyperparameter – which is what
[`distrib_cross_y`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.html)
exists for.

## Usage

``` r
distrib_penalty(d, map = NULL, n_coef = NULL, kinks = numeric(0))
```

## Arguments

- d:

  A univariate continuous distributions7 object; typically a `fixed()`
  wrapper holding the location at zero.

- map:

  The matrix \\D\\, or `NULL` (default) for the identity.

- n_coef:

  The number of coefficients; required when `map` is `NULL`, ignored
  otherwise.

- kinks:

  The points where the parent's log-density is not differentiable in its
  argument, declared by the caller because a distribution object does
  not carry a response kink set; `0` for a Laplace at zero.

## Value

An object of class `DistribPenalty`.

## Details

Ridge, lasso and the heavy-tailed prior are this construction at a
[`fixed`](https://statmodels7.github.io/distributions7/reference/fixed.html)
Gaussian, Laplace and Student t; see
[`ridge_penalty`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md).
The normalizing constant comes with the density and is kept, so the
value is exactly the negative log prior density and a free scale or a
free \\\nu\\ is estimable.

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
```
