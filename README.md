# penalties7 <img src="man/figures/logo.png" align="right" height="139" alt="" />

Penalties for regularized and Bayesian regression as S7 objects. A penalty
is rho(D beta; theta) -- a linear map, a scalar function, hyperparameters --
and the object answers everything an estimation routine asks: the value
with its normalizing constant, the exact derivatives in the coefficients
and in the hyperparameters, the mixed block, the kink set, and, for
quadratic penalties, the rank, the null basis and the log
pseudo-determinant a marginal criterion consumes.

Three branches: `quadratic_penalty()`; `distrib_penalty()`, a univariate
[distributions7](https://statmodels7.github.io/distributions7/) log-density
applied coordinatewise (`ridge_penalty()`, `lasso_penalty()`,
`heavy_penalty()` are its named instances); and `scad_penalty()` /
`mcp_penalty()`, defined by their derivative and improper by construction.
Part of the [statmodels7](https://statmodels7.github.io) toolkit.

## Installation

``` r
# install.packages("pak")
pak::pak("statmodels7/penalties7")
```
