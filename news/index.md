# Changelog

## penalties7 0.1.0

- First release: penalties as S7 objects, rho(D beta; theta), with the
  value, the exact derivatives in the coefficients and the
  hyperparameters, the mixed block, and the declared kink set. Three
  branches:
  [`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
  (rank, null basis and log pseudo-determinant fixed at construction,
  the pieces a marginal criterion consumes),
  [`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md)
  (a univariate distributions7 log-density applied coordinatewise, with
  [`ridge_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md),
  [`lasso_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md)
  and
  [`heavy_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md)
  as named instances), and
  [`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)/[`mcp_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
  (defined by their derivative, improper by construction). Normalizing
  constants are kept throughout, so a proper penalty is exactly the
  negative log-density of its prior.
  [`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
  validates every closed form against a route sharing no code with it.
