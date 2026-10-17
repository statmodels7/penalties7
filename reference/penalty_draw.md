# Drawing From a Penalty Read as a Prior

Returns a draw from the distribution whose negative log-density is the
penalty, at the given hyperparameters, with `NA` in the coordinates that
the distribution does not determine.

## Usage

``` r
penalty_draw(pen, theta, ...)
```

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- theta:

  A named list of hyperparameter values, or a named numeric vector
  carrying the same, as for the other generics of the package.

- ...:

  Passed to methods. The shipped methods do not read it.

## Value

A numeric vector of length `pen@n_coef`, `NA` in the coordinates that
the prior does not determine and possibly `NA` throughout.

## Details

A proper penalty is a negative log-prior, and
[`penalty_value()`](https://statmodels7.github.io/penalties7/reference/penalty_value.md)
includes its normalizing constant for that reason. Simulating from a
model that contains a random effect, a ridge or a lasso therefore
amounts to drawing those coefficients from this distribution.

## Branches that draw

A separable penalty draws from its own parent family, so a Gaussian
prior gives Gaussian effects, a Laplace prior gives Laplace effects and
a Student t prior gives heavy-tailed effects. A univariate parent is
drawn once per coordinate, and a multivariate parent once per block of
coordinates.

A quadratic, additive or structured penalty is a Gaussian prior with
precision \\S(\theta)\\, and two cases are drawn directly. Where \\S\\
is diagonal each coordinate is drawn on its own with standard deviation
\\S\_{jj}^{-1/2}\\, which covers a ridge and a smooth in a basis whose
penalty is diagonal (such as the Demmler-Reinsch basis); a zero on that
diagonal is a direction that the prior does not determine, and its entry
is `NA`. Where \\S\\ is not diagonal but has full rank, one Cholesky
factor \\S = R^\top R\\ gives the draw, since \\\beta = R^{-1}z\\ has
covariance \\S^{-1}\\.

A rank-deficient \\S\\ that is not diagonal, such as the penalty of an
anisotropic tensor smooth, gives a vector of `NA`. The prior is flat
along the null space of \\S\\, and a draw restricted to the range would
require a basis of that range, a decomposition that these branches do
not compute.

SCAD and MCP are improper and are not densities, so their draws are `NA`
everywhere. A branch that is not quadratic and does not register a
method of its own also gives `NA` everywhere.

## The map

A separable penalty is a prior on \\D\beta\\, so a draw of the
coefficients inverts the map. Without a map the draw is returned as it
is. Under a diagonal map, such as the one that standardization produces,
each coordinate is divided by its diagonal entry, and a coordinate whose
entry is zero is not determined by the prior and is `NA`. Under any
other map, including a square invertible one, the draw is not inverted
and every coordinate is `NA`. A quadratic penalty does not need this
step, because
[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
returns the precision in the coefficients themselves.

## Undetermined coordinates

The coordinates that the prior does not determine are `NA` rather than
zero. A caller that draws a whole model has its own rule for a
coefficient that is not covered by any prior; an intercept is the
ordinary case. A zero would be indistinguishable from a prior that
concentrates at zero, so it would not show which coordinates remain to
be filled.

## Methods

The method on the base class
[`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
draws from a quadratic penalty (for which
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
is `TRUE`) read as a Gaussian prior, as described above, and returns
`NA` for any other branch. A separable penalty registers
[`penalty_draw.DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_draw.DistribPenalty.md).

## See also

[`penalty_value()`](https://statmodels7.github.io/penalties7/reference/penalty_value.md)
for the negative log-density from which the draw is made,
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
for whether the penalty is a proper prior,
[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
for the precision of the Gaussian branches.

## Examples

``` r
set.seed(1)

# A Gaussian prior on ten effects gives ten Gaussian draws.
ridge <- ridge_penalty(n_coef = 10)
round(penalty_draw(ridge, list(lambda = 4)), 3)
#>  [1] -0.313  0.092 -0.418  0.798  0.165 -0.410  0.244  0.369  0.288 -0.153

# The scale is that of the prior, so a larger lambda shrinks the draw.
c(loose = stats::sd(penalty_draw(ridge, list(lambda = 0.25))),
  tight = stats::sd(penalty_draw(ridge, list(lambda = 100))))
#>      loose      tight 
#> 2.13902967 0.09556076 

# The lasso is a Laplace prior and gives Laplace draws.
round(penalty_draw(lasso_penalty(n_coef = 6), list(lambda = 2)), 3)
#> [1]  0.874 -0.266 -0.043 -0.204  0.180 -0.331

# The unpenalized direction of a smooth is NA.
P <- diag(c(0, rep(1, 4)))
penalty_draw(quadratic_penalty(P), list(lambda = 1))
#> [1]          NA -0.05380504 -1.37705956 -0.41499456 -0.39428995

# A diagonal map with a zero entry leaves that coordinate NA.
penalty_draw(lasso_penalty(map = Matrix::Diagonal(x = c(1, 0, 2))),
             list(lambda = 1))
#> [1] -0.04845261          NA  0.65222622

# SCAD is not a density, so every coordinate is NA.
penalty_draw(scad_penalty(n_coef = 3), list(lambda = 1, a = 3.7))
#> [1] NA NA NA
```
