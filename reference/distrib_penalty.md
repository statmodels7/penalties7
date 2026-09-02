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
  wrapper holding the location at zero. A discrete or otherwise
  non-continuous object is rejected, naming its `dimension`.

- map:

  The matrix \\D\\, with one column per coefficient, or `NULL` (the
  default) for the identity. A Matrix object keeps its own storage, a
  diagonal map being what standardization is.

- n_coef:

  The number of coefficients. Required when `map` is `NULL` and ignored
  otherwise, `ncol(map)` being the count then.

- kinks:

  The points where the parent's log-density is not differentiable in its
  argument. `NULL`, the default, derives them with
  [`distrib_kinks()`](https://statmodels7.github.io/penalties7/reference/distrib_kinks.md);
  pass a numeric vector to say so directly, or `numeric(0)` to declare
  there are none. A multivariate parent is given `numeric(0)` without
  asking, a kink being a point of a scalar argument.

## Value

A
[`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
object whose hyperparameters, bounds and links are the parent
distribution's.

## The hyperparameters are the distribution's

`params`, `params_bounds` and `link_params` are read off the
distribution rather than restated, and every derivative is the
distribution's own, reassembled through the map:

\$\$\frac{\partial\rho}{\partial\beta} = -D'\ell^{(y)}, \qquad
\frac{\partial^2\rho}{\partial\beta^2} =
-D'\\\mathrm{diag}(\ell^{(yy)})\\D, \qquad
\frac{\partial^2\rho}{\partial\beta\\\partial\theta_k} =
-D'\ell^{(y\theta_k)},\$\$

with the middle matrix block diagonal instead of diagonal when the
parent is multivariate. The hyperparameter blocks are the parent's score
and Hessian summed over the blocks, and the mixed block is
[`distributions7::distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.html),
which exists for this.

## The constant is kept

The normalizing constant comes with the density and is not dropped, so
the value is exactly the negative log prior density. A free scale or a
free \\\nu\\ is then estimable: with the constant dropped, a prior scale
could be sent to infinity for nothing.

## A multivariate parent

Centering is the caller's, typically through
[`distributions7::fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.html)
at a zero mean, and the matrix parameter carries the dependence within a
block. The number of values the penalty is read at must divide into
whole blocks, and the constructor rejects a count that does not, naming
both numbers.

A blockwise penalty has **no proximal operator**: that operator acts one
coordinate at a time and the coordinates of a block do not separate.
[`has_prox()`](https://statmodels7.github.io/penalties7/reference/has_prox.md)
answers `FALSE` and
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
rejects.

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
# A Gaussian prior at zero: the separable twin of the ridge, written by
# its standard deviation rather than by a precision.
d <- distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0)
pen <- distrib_penalty(d, n_coef = 3)
pen@params
#> [1] "sigma"
penalty_value(pen, c(1, 0, -1), list(sigma = 2))
#> [1] 5.086257

# Which is exactly the negative log-density of that prior.
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

# And it has no proximal operator, the coordinates of a block being tied.
has_prox(pen2)
#> [1] FALSE

# A count that does not divide into blocks is rejected.
try(distrib_penalty(mv, n_coef = 5))
#> Error : 'fixed multivariate gaussian [2d, sigma=log_cholesky] [mu1=0,mu2=0]' is 2-variate, so the 5 values it is read at must divide into
#>   whole blocks; 5 does not.
```
