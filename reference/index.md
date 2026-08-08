# Package index

## The class and the generics

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

## The quadratic branch

- [`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
  : Construct a Quadratic Penalty
- [`QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/QuadraticPenalty.md)
  : S7 Class for the Quadratic Penalty

## The separable branch

- [`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md)
  : Construct a Separable Penalty From a Distribution
- [`ridge_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md)
  [`lasso_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md)
  [`heavy_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md)
  : Named Separable Penalties
- [`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
  : S7 Class for the Separable Penalty

## The structured branch

- [`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md)
  : Construct a Structured Quadratic Penalty
- [`StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/StructuredPenalty.md)
  : S7 Class for the Structured Quadratic Penalty

## SCAD and MCP

- [`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
  [`mcp_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
  : Construct the SCAD and MCP Penalties
- [`ScadPenalty()`](https://statmodels7.github.io/penalties7/reference/ScadPenalty.md)
  [`McpPenalty()`](https://statmodels7.github.io/penalties7/reference/ScadPenalty.md)
  : S7 Classes for the Derivative-Defined Penalties

## The proximal operator

- [`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
  : Proximal Operator of a Penalty
- [`has_prox()`](https://statmodels7.github.io/penalties7/reference/has_prox.md)
  : Does a Penalty Supply a Proximal Operator?

## Validation

- [`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
  : Check a Penalty Numerically
