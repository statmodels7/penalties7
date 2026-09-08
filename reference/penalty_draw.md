# Drawing From a Penalty Read as a Prior

A vector drawn from the distribution whose negative log-density the
penalty is, at the hyperparameters given, with `NA` in the coordinates
that distribution does not determine.

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
  carrying the same, as every other generic here takes them.

- ...:

  Passed to methods. No shipped method reads it.

## Value

A numeric vector of length `pen@n_coef`, `NA` in the coordinates the
prior does not determine and possibly `NA` throughout.

## Details

A penalty is a negative log-prior, which is why
[`penalty_value()`](https://statmodels7.github.io/penalties7/reference/penalty_value.md)
keeps the normalizing constant, and a prior is something one draws from.
Simulating from a model that carries a random effect, a ridge or a lasso
means drawing those coefficients from exactly this distribution instead
of from a normal chosen by whoever wrote the simulation.

## Which branches answer

A separable penalty draws coordinatewise from its own parent family, so
a Gaussian prior gives Gaussian effects, a Laplace prior gives Laplace
ones and a Student t prior gives heavy-tailed ones. That is the case a
random effect is, and the one this exists for.

A quadratic, additive or structured penalty is a Gaussian prior with
precision \\S(\theta)\\, and two of its shapes are cheap. Where \\S\\ is
diagonal each coordinate is drawn on its own at a standard deviation
\\S\_{jj}^{-1/2}\\, which covers a ridge and a Demmler-Reinsch smooth,
whose penalty is \\\mathrm{diag}(0, 1, \dots, 1)\\; a zero on that
diagonal is a direction the prior says nothing about and comes back
`NA`. Where \\S\\ is not diagonal but has full rank, one Cholesky factor
gives the draw, \\\beta = R^{-1}z\\ having covariance \\S^{-1}\\.

A deficient non-diagonal \\S\\, which an anisotropic tensor smooth has,
comes back all `NA`. The prior is flat along its null space and drawing
on the range alone needs a basis of that range, which is a decomposition
of a matrix these branches are built not to form.

SCAD and MCP are improper by construction, densities of nothing, and
answer `NA` everywhere. So does any penalty whose branch says nothing
here.

## The map

The penalty is a prior on \\D\beta\\, so a draw of the coefficients
inverts the map. For a separable penalty the identity and a diagonal
map, which is what standardization is, invert coordinatewise, and
anything else comes back `NA`: under a general map the prior leaves the
coefficients underdetermined, which is the same fact that makes
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
reject one. A quadratic penalty needs no such care,
[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
returning the precision in the coefficients themselves.

## Why the unanswered coordinates are `NA` rather than zero

A caller drawing a whole model has a rule of its own for a coefficient
no prior covers, an intercept being the ordinary case. Returning zero
would be indistinguishable from a prior that genuinely concentrates
there, and the caller could not tell which coordinates it still has to
fill.

## See also

[`penalty_value()`](https://statmodels7.github.io/penalties7/reference/penalty_value.md)
for the log-density this inverts,
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
for whether there is a distribution to draw from at all,
[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
for the precision the Gaussian branches read.

## Examples

``` r
set.seed(1)

# A Gaussian prior on ten effects: ten Gaussian effects.
ridge <- ridge_penalty(n_coef = 10)
round(penalty_draw(ridge, list(lambda = 4)), 3)
#>  [1] -0.313  0.092 -0.418  0.798  0.165 -0.410  0.244  0.369  0.288 -0.153

# The scale is the prior's own, so a larger lambda shrinks the draw.
c(loose = stats::sd(penalty_draw(ridge, list(lambda = 0.25))),
  tight = stats::sd(penalty_draw(ridge, list(lambda = 100))))
#>      loose      tight 
#> 2.13902967 0.09556076 

# A Laplace prior gives Laplace effects, which is what a lasso says.
round(penalty_draw(lasso_penalty(n_coef = 6), list(lambda = 2)), 3)
#> [1]  0.874 -0.266 -0.043 -0.204  0.180 -0.331

# A smooth leaves its null space undetermined, and says so.
P <- diag(c(0, rep(1, 4)))
penalty_draw(quadratic_penalty(P), list(lambda = 1))
#> [1]          NA -0.05380504 -1.37705956 -0.41499456 -0.39428995

# SCAD is the density of nothing.
penalty_draw(scad_penalty(n_coef = 3), list(lambda = 1, a = 3.7))
#> [1] NA NA NA
```
