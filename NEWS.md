# penalties7 0.3.0

* penalty_prox(): the proximal operator, closed form for the quadratic
  and structured branches (one linear solve), for the Gaussian and
  Laplace instances, and for SCAD and MCP over their piecewise
  regions; any other separable penalty is solved coordinatewise from
  its response derivative. has_prox() asks before calling.

# penalties7 0.2.0

* The structured quadratic prior: `structured_penalty()` takes a
  parameters7 matrix_parameter as the PRECISION, so the hyperparameters
  reach every entry of the matrix. The free vector is unconstrained by
  construction, so every link is the identity -- the flattening convention
  of the multivariate families. Every derivative comes from the
  structure's own contract (param_d1/param_d2, the logdet derivatives),
  and the marginal pieces answer through it, so is_quadratic() is TRUE.
  At a zero free vector the log-Cholesky prior IS the plain ridge, pinned
  at machine precision. check_penalty() gains the structured logpdet
  check, its lambda-slope identity now gated on the branch that has a
  lambda.

# penalties7 0.1.0

* First release: penalties as S7 objects, rho(D beta; theta), with the
  value, the exact derivatives in the coefficients and the
  hyperparameters, the mixed block, and the declared kink set. Three
  branches: `quadratic_penalty()` (rank, null basis and log
  pseudo-determinant fixed at construction, the pieces a marginal
  criterion consumes), `distrib_penalty()` (a univariate distributions7
  log-density applied coordinatewise, with `ridge_penalty()`,
  `lasso_penalty()` and `heavy_penalty()` as named instances), and
  `scad_penalty()`/`mcp_penalty()` (defined by their derivative, improper
  by construction). Normalizing constants are kept throughout, so a
  proper penalty is exactly the negative log-density of its prior.
  `check_penalty()` validates every closed form against a route sharing
  no code with it.
