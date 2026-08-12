# Changelog

## penalties7 0.9.0

- [`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
  describes the scalar proximal operator of a separable penalty as an
  odd piecewise linear table, so a compiled loop can apply it without
  knowing which family it came from. Every closed form the package
  carries has that shape: the soft threshold is two pieces, the elastic
  net two, MCP three, SCAD four, and a Gaussian prior one. The step is a
  vector, one per coefficient, because in a coordinate descent the step
  of coordinate j is 1/sum(w x_j^2) and does not move while the working
  weights are held.

  The operator is applied once per coordinate per sweep at a point that
  moves every time, so a compiled loop calling back for it would spend
  the gain on the calls. Passing the numbers instead keeps the
  mathematics in the penalty and leaves the kernel naming no family.

  [`prox_apply()`](https://statmodels7.github.io/penalties7/reference/prox_apply.md)
  evaluates a table in R, and the table is pinned against
  [`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
  itself across every breakpoint and at the breakpoints exactly, at
  three step lengths and five families. A penalty with no such
  description – a quadratic under a general matrix, an operator that is
  a root rather than a formula, a parent not centred where the quadratic
  pull is, a step past the convex region of SCAD or MCP – returns
  `NULL`.

## penalties7 0.8.0

- [`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md)
  derives its kink set from the parent instead of defaulting to none. A
  penalty built by hand from a non-smooth family –
  `distrib_penalty(fixed(laplace_distrib(), mu = 0))`, which is the
  lasso – declared itself differentiable everywhere, so a model layer
  reading
  [`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
  put it in the scheme for the opposite property. The shipped instances
  passed `kinks` explicitly and were never affected.

  [`distrib_kinks()`](https://statmodels7.github.io/penalties7/reference/distrib_kinks.md)
  takes the candidates from `params_smooth` crossed with what `fixed()`
  holds, a location that is not smooth being a kink in the argument at
  the value it is held at, and then measures each one by comparing the
  one-sided derivatives of the log-density across it. Inferring alone
  would put a kink on any family whose non-smooth parameter is not a
  location; measuring alone would need somewhere to look. Nothing is
  taken from a parameter that is free, its value being whatever the
  hyperparameters say at the time.

  Passing `kinks` still overrides, including `numeric(0)` to declare
  there are none.

## penalties7 0.7.1

- The hyperparameters may be given as a named numeric vector as well as
  the documented list, the alignment converting one to the other. The
  branches had split on how they read `theta`: `[[` accepts both shapes
  and `$` accepts only the list, so a caller passing a vector reached
  the quadratic and separable branches and stopped inside `scad()` and
  `mcp()` with “\$ operator is invalid for atomic vectors”, three frames
  down and naming neither the argument nor the penalty.

## penalties7 0.7.0

- [`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md)
  and
  [`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md)
  read the parent’s `distrib_grad_y_hess()` and `distrib_hess_y_hess()`
  for the separable branch instead of differencing its first-order
  components. With a gaussian parent – every ridge, every random effect
  – the branch is now exact: measured, the second derivative of
  `I/sigma^2` comes back as `6I/sigma^4` to 1e-13 where the difference
  gave 1e-10.

## penalties7 0.6.0

- [`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md),
  [`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md)
  and
  [`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md)
  are what a marginal criterion asks of a penalty: the hyperparameter
  derivatives of the coefficient Hessian and of the mixed block. Every
  branch answers – the quadratic and additive ones from their
  components, the structured one from the matrix parameter’s `param_d1`
  and `param_d2`, the separable one from the parent’s response surface –
  so a penalty is estimable by REML or ML whatever its shape, and one
  with a kink rejects by name.
  [`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md)
  reports whether the third derivative in the coefficients is zero,
  asking the parent whether its response curvature depends on the
  response rather than recognizing a family by name.

## penalties7 0.5.0

- elasticnet_penalty(): the elastic net as the same construction as
  ridge and lasso, a separable penalty over a fixed() family – here
  distributions7::enet_distrib(), the product of the Laplace and the
  Gaussian at zero, normalized. Its hyperparameters are the overall rate
  lambda and the mixing weight alpha, and the normalizing constant
  depends on both, which is what makes them estimable by a marginal
  criterion.
- penalty_prox() gains the elastic-net closed form, the soft threshold
  of the Laplace part followed by the shrinkage of the Gaussian one.

## penalties7 0.4.0

- additive_penalty(): a sum of quadratic penalties with a smoothing
  parameter on each, which is what an anisotropic tensor smooth needs.
  The rank is fixed at construction from the components stacked and
  normalized, since the null space of the sum is the intersection of
  theirs and does not move with the parameters – a count taken from the
  assembled matrix does.

## penalties7 0.3.0

- penalty_prox(): the proximal operator, closed form for the quadratic
  and structured branches (one linear solve), for the Gaussian and
  Laplace instances, and for SCAD and MCP over their piecewise regions;
  any other separable penalty is solved coordinatewise from its response
  derivative. has_prox() asks before calling.

## penalties7 0.2.0

- The structured quadratic prior:
  [`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md)
  takes a parameters7 matrix_parameter as the PRECISION, so the
  hyperparameters reach every entry of the matrix. The free vector is
  unconstrained by construction, so every link is the identity – the
  flattening convention of the multivariate families. Every derivative
  comes from the structure’s own contract (param_d1/param_d2, the logdet
  derivatives), and the marginal pieces answer through it, so
  is_quadratic() is TRUE. At a zero free vector the log-Cholesky prior
  IS the plain ridge, pinned at machine precision. check_penalty() gains
  the structured logpdet check, its lambda-slope identity now gated on
  the branch that has a lambda.

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
