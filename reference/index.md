# Package index

## The class and the generics

What every penalty answers: its value with the normalizing constant, its
derivatives in the coefficients and in the hyperparameters, the mixed
block, and the kink set a non-smooth method needs.

- [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  : S7 Base Class for Penalties
- [`penalty_value()`](https://statmodels7.github.io/penalties7/reference/penalty_value.md)
  : Value of a Penalty
- [`penalty_gradient()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
  [`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
  : Coefficient Derivatives of a Penalty
- [`penalty_grad_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
  [`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
  [`penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
  : Hyperparameter Derivatives of a Penalty
- [`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
  : The Non-Differentiable Points of a Penalty
- [`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
  : Is a Penalty a Proper Prior?
- [`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
  : Is a Penalty Quadratic?
- [`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
  [`penalty_rank()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
  [`penalty_null_basis()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
  [`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
  : The Pieces a Marginal Criterion Consumes

## What a marginal criterion asks

A Laplace approximation to a marginal likelihood differentiates the
determinant of the penalized information, so it needs the hyperparameter
derivatives of the coefficient Hessian and of the mixed block. A penalty
that answers these is estimable by REML or ML whatever its shape.

- [`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
  : The Derivative of the Coefficient Hessian in the Hyperparameters
- [`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md)
  : The Second Derivative of the Coefficient Hessian in the
  Hyperparameters
- [`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md)
  : The Derivative of the Mixed Block in the Hyperparameters
- [`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md)
  : Is a Penalty Quadratic in the Coefficients?

## The quadratic branch

A fixed matrix and one smoothing parameter, so the penalty is the
negative log-density of an improper Gaussian prior. The rank, the null
basis and the log pseudo-determinant are fixed at one eigendecomposition
when the object is built.

- [`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
  : Construct a Quadratic Penalty
- [`QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/QuadraticPenalty.md)
  : S7 Class for the Quadratic Penalty

## The additive branch

A sum of quadratics with a smoothing parameter on each, which is what an
anisotropic tensor smooth needs. The rank comes from the components
stacked and normalized, not from the assembled matrix, whose eigenvalue
count falls as the parameters spread apart.

- [`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
  : Construct a Sum of Quadratic Penalties
- [`AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/AdditivePenalty.md)
  : S7 Class for a Sum of Quadratic Penalties

## The separable branch

A univariate distributions7 log-density applied coordinatewise. Ridge is
a Gaussian at zero, the lasso a Laplace, the elastic net the product of
the two, and the heavy-tailed prior a Student t.

- [`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md)
  : Construct a Separable Penalty From a Distribution
- [`ridge_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md)
  [`lasso_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md)
  [`elasticnet_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md)
  [`heavy_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md)
  : Named Separable Penalties
- [`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
  : S7 Class for the Separable Penalty

## The structured branch

A parameters7 matrix parameter used as the precision, so the
hyperparameters are that structure’s free values and every derivative
comes from its derivative arrays.

- [`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md)
  : Construct a Structured Quadratic Penalty
- [`StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/StructuredPenalty.md)
  : S7 Class for the Structured Quadratic Penalty

## SCAD and MCP

Defined by their derivative rather than by a density, hence improper by
construction: they have no normalizing constant and cannot come from a
distribution.

- [`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
  [`mcp_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
  : Construct the SCAD and MCP Penalties
- [`ScadPenalty()`](https://statmodels7.github.io/penalties7/reference/ScadPenalty.md)
  [`McpPenalty()`](https://statmodels7.github.io/penalties7/reference/ScadPenalty.md)
  : S7 Classes for the Derivative-Defined Penalties

## The proximal operator

One step of a proximal gradient method, in closed form where one exists
and as a linear solve or a coordinatewise root otherwise.

- [`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
  : Proximal Operator of a Penalty
- [`has_prox()`](https://statmodels7.github.io/penalties7/reference/has_prox.md)
  : Does a Penalty Supply a Proximal Operator?

## Validation

The numerical checks a penalty must pass, meant above all for one
written outside the package.

- [`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
  : Check a Penalty Numerically
