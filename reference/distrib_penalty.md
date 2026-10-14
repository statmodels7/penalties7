# Construct a Separable Penalty From a Distribution

Builds the penalty \\\rho(D\beta;\theta) = -\sum_i \log f(b_i;\theta)\\,
where \\f\\ is a distributions7 density and the \\b_i\\ are the
successive blocks of \\D\beta\\. A univariate parent gives blocks of one
coordinate, which is the separable penalty; a \\p\\-variate parent gives
blocks of \\p\\, a prior under which the coordinates of one block depend
on each other while the blocks stay independent.

## Usage

``` r
distrib_penalty(d, map = NULL, n_coef = NULL, kinks = NULL)
```

## Arguments

- d:

  A continuous distributions7 object, univariate or multivariate;
  typically a
  [`distributions7::fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.html)
  wrapper holding the location at zero. A univariate object must inherit
  from
  [distributions7::continuous_distrib](https://statmodels7.github.io/distributions7/reference/continuous_distrib.html);
  otherwise the constructor signals an error that names the
  distribution. A multivariate object is rejected in the same way when
  it is, or wraps, the distributions7 multinomial; a multivariate family
  written outside distributions7 is not checked. An object whose
  `dimension` is neither `"univariate"` nor `"multivariate"` is
  rejected.

- map:

  The matrix \\D\\, with one column per coefficient, or `NULL` (the
  default) for the identity. A Matrix object keeps its storage, so a
  diagonal or sparse map is not densified, and a square map whose
  off-diagonal entries are all zero is stored as a
  [`Matrix::Diagonal()`](https://rdrr.io/pkg/Matrix/man/Diagonal.html).

- n_coef:

  The number of coefficients. It is required when `map` is `NULL` and
  ignored otherwise, the count being `ncol(map)`.

- kinks:

  The points where the parent's log-density is not differentiable in its
  argument. With `NULL`, the default, they are derived by
  [`distrib_kinks()`](https://statmodels7.github.io/penalties7/reference/distrib_kinks.md)
  for a univariate parent and set to `numeric(0)` for a multivariate
  one, a kink being a point of a scalar argument. A numeric vector is
  used as given, and `numeric(0)` declares that there are none.

## Value

A
[`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
object whose hyperparameters, bounds, links and `params_smooth` are
those of the parent distribution.

## Hyperparameters and derivatives

The properties `params`, `params_bounds`, `link_params` and
`params_smooth` are copied from the distribution, and every derivative
is the distribution's own, carried through the map:

\$\$\frac{\partial\rho}{\partial\beta} = -D'\ell^{(y)}, \qquad
\frac{\partial^2\rho}{\partial\beta^2} =
-D'\\\mathrm{diag}(\ell^{(yy)})\\D, \qquad
\frac{\partial^2\rho}{\partial\beta\\\partial\theta_k} =
-D'\ell^{(y\theta_k)},\$\$

where the middle matrix is block diagonal instead of diagonal when the
parent is multivariate. The hyperparameter derivatives are the parent's
score and Hessian summed over the blocks, and the mixed derivative comes
from
[`distributions7::distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.html).

## Normalizing constant and Jacobian

The normalizing constant of the density is kept, so the value is the
negative log prior density. This makes a free scale or a free \\\nu\\
estimable by a marginal criterion. With the constant dropped, the value
of a Gaussian prior would decrease toward zero as its scale grows, and a
criterion would send the scale to infinity.

Under a square nonsingular map the value also carries the Jacobian of
\\\beta \mapsto D\beta\\: it is \\-\sum_i \log f(b_i;\theta) -
\log\lvert\det D\rvert\\, the negative log-density of \\\beta\\ instead
of that of \\D\beta\\. Under a diagonal map, which is the map that
standardization builds, the term is \\\sum_j \log\lvert d_j\rvert\\. The
term is constant in \\\beta\\ and in \\\theta\\, so it changes no
derivative, and it makes a marginal criterion invariant to the units of
a covariate. The determinant is computed once, when the penalty is
built, and stored in the `log_jacobian` property; for a diagonal map it
is recomputed from the entries each time the value is computed. Under a
map that is not square, or is singular, \\\exp(-\rho)\\ does not
integrate to one over \\\beta\\,
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
returns `FALSE`, and the value is \\-\sum_i \log f(b_i;\theta)\\ with no
further constant.

A square base matrix or Matrix object whose off-diagonal entries are all
exactly zero is stored as a
[`Matrix::Diagonal()`](https://rdrr.io/pkg/Matrix/man/Diagonal.html)
map. The proximal operator, the draws and the Jacobian recognize a
diagonal map by its class, so this conversion makes their result depend
on the entries of the map and not on the class in which it was supplied.

## A multivariate parent

The caller centers the parent, typically through
[`distributions7::fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.html)
at a zero mean, and the matrix parameter carries the dependence within a
block. The number of values at which the penalty is read must divide
into whole blocks; otherwise the constructor signals an error that names
the block width and the count.

A blockwise penalty has **no proximal operator**, because that operator
acts one coordinate at a time and the coordinates of a block do not
separate.
[`has_prox()`](https://statmodels7.github.io/penalties7/reference/has_prox.md)
returns `FALSE` and
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
signals an error.

## See also

[`ridge_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md),
[`lasso_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md),
[`elasticnet_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md)
and
[`heavy_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md)
for the named instances,
[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
and
[`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md)
for the quadratic branches,
[`distrib_kinks()`](https://statmodels7.github.io/penalties7/reference/distrib_kinks.md)
for how the kinks are found,
[`distributions7::fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.html)
for the wrapper that centers a parent.

## Examples

``` r
# A Gaussian prior at zero: the separable form of the ridge, written by
# its standard deviation instead of a precision.
d <- distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0)
pen <- distrib_penalty(d, n_coef = 3)
pen@params
#> [1] "sigma"
penalty_value(pen, c(1, 0, -1), list(sigma = 2))
#> [1] 5.086257

# The same number from stats::dnorm().
-sum(stats::dnorm(c(1, 0, -1), sd = 2, log = TRUE))
#> [1] 5.086257

# A correlated prior over two coefficients per block, three blocks. The
# hyperparameters are the matrix parameter's free values.
mv <- distributions7::fixed(distributions7::mvgaussian1_distrib(2),
                            mu1 = 0, mu2 = 0)
pen2 <- distrib_penalty(mv, n_coef = 6)
pen2@block
#> [1] 2
pen2@params
#> [1] "sigma_log_L1" "sigma_log_L2" "sigma_L2.1"  
penalty_value(pen2, c(1, 0, -1, 0.5, 0.2, -0.3),
              list(sigma_log_L1 = 0, sigma_log_L2 = 0, sigma_L2.1 = 0.4))
#> [1] 7.090831

# The blockwise penalty has no proximal operator.
has_prox(pen2)
#> [1] FALSE

# A count that does not divide into blocks is rejected.
try(distrib_penalty(mv, n_coef = 5))
#> Error : 'fixed multivariate gaussian [2d, sigma=log_cholesky] [mu1=0,mu2=0]' is 2-variate, so the 5 values it is read at must divide into
#>   whole blocks; 5 does not.

# A square nonsingular map adds -log|det D| to the value, and a diagonal
# base matrix is stored as a Matrix::Diagonal.
D <- matrix(c(2, 1, 0, 3), 2)
pen3 <- distrib_penalty(d, map = D)
penalty_value(pen3, c(0.4, -0.7), list(sigma = 1))
#> [1] 1.811118
-sum(stats::dnorm(D %*% c(0.4, -0.7), log = TRUE)) - log(abs(det(D)))
#> [1] 1.811118
methods::is(distrib_penalty(d, map = diag(c(2, 3)))@map, "diagonalMatrix")
#> [1] TRUE
```
