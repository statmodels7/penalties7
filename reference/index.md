# Package index

## The class and the generics

The base class and the generics of every penalty: the value with its
normalizing constant, the derivatives in the coefficients and in the
hyperparameters, the mixed block, the points where the penalty is not
smooth, and the predicates and quantities that a fitting layer reads.

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
  : Points Where a Penalty Is Not Smooth
- [`penalty_readable()`](https://statmodels7.github.io/penalties7/reference/penalty_readable.md)
  : Interpretable Quantities of a Penalty's Hyperparameters
- [`penalty_draw()`](https://statmodels7.github.io/penalties7/reference/penalty_draw.md)
  : Drawing From a Penalty Read as a Prior
- [`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
  : Is a Penalty a Proper Prior?
- [`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
  : Is a Penalty a Quadratic Form?
- [`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
  [`penalty_rank()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
  [`penalty_null_basis()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
  [`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
  : Quantities of a Quadratic Penalty for a Marginal Criterion

## Derivatives for a marginal criterion

A Laplace approximation to a marginal likelihood differentiates the
determinant of the penalized information, so it needs the hyperparameter
derivatives of the coefficient Hessian and of the mixed block. The
hyperparameters of a penalty that supplies them can be estimated by REML
or ML. SCAD, MCP and a separable penalty with a kink do not supply them.

- [`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
  : The Derivative of the Coefficient Hessian in the Hyperparameters
- [`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md)
  : The Second Derivative of the Coefficient Hessian in the
  Hyperparameters
- [`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md)
  : The Derivative of the Mixed Block in the Hyperparameters
- [`penalty_dhessian_beta()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian_beta.md)
  : The Derivative of the Coefficient Hessian in the Coefficients
- [`penalty_d2hessian_beta()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.md)
  [`penalty_dhessian_beta_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.md)
  : Second Derivatives of the Coefficient Hessian in the Coefficients
  and the Hyperparameters
- [`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md)
  : Is a Penalty Quadratic in the Coefficients?

## The quadratic branch

A fixed matrix P, an optional map D and one smoothing parameter lambda,
so the penalty is the negative log-density of a Gaussian prior with
precision lambda D’PD, improper when D’PD is singular. The rank, the
null basis and the log pseudo-determinant are computed from D’PD when
the object is built.

- [`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
  : Construct a Quadratic Penalty
- [`QuadraticPenalty()`](https://statmodels7.github.io/penalties7/reference/QuadraticPenalty.md)
  : S7 Class for the Quadratic Penalty

## The additive branch

A sum of quadratic penalties with a smoothing parameter on each, which
is the penalty of an anisotropic tensor-product smooth. The rank is
computed at construction from the sum of the components, each divided by
its largest absolute entry, because the null space of the sum does not
depend on the smoothing parameters.

- [`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
  : Construct a Sum of Quadratic Penalties
- [`AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/AdditivePenalty.md)
  : S7 Class for a Sum of Quadratic Penalties

## The separable branch

The negative log-density of a continuous distributions7 family at the
coordinates of the mapped coefficients, one coordinate at a time for a
univariate family and one block at a time for a multivariate one. The
lasso is a Laplace at zero, the elastic net the product of a Laplace and
a Gaussian, and the heavy-tailed prior a Student t. The ridge,
documented with them, is the Gaussian at zero written by its precision
and is built on the quadratic branch.

- [`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md)
  : Construct a Separable Penalty From a Distribution
- [`distrib_kinks()`](https://statmodels7.github.io/penalties7/reference/distrib_kinks.md)
  : Where the Parent's Log-Density Has a Kink
- [`ridge_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md)
  [`lasso_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md)
  [`elasticnet_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md)
  [`heavy_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md)
  : The Named Penalties
- [`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
  : S7 Class for the Separable Penalty

## The structured branch

A parameters7 matrix parameter used as the precision, so the
hyperparameters are the free values of that structure and every
derivative comes from its derivative arrays.

- [`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md)
  : Construct a Structured Quadratic Penalty
- [`StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/StructuredPenalty.md)
  : S7 Class for the Structured Quadratic Penalty

## SCAD and MCP

Defined by their derivative rather than by a density. Their value is
bounded, so exp(-value) is not integrable: they are improper, carry no
normalizing constant and cannot come from a distribution.

- [`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
  [`mcp_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
  : Construct the SCAD and MCP Penalties
- [`ScadPenalty()`](https://statmodels7.github.io/penalties7/reference/ScadPenalty.md)
  [`McpPenalty()`](https://statmodels7.github.io/penalties7/reference/ScadPenalty.md)
  : S7 Classes for the SCAD and MCP Penalties

## The proximal operator

The operator that a proximal gradient method applies after each gradient
step: a linear solve for the quadratic and structured branches, a closed
form for the Gaussian, Laplace, elastic-net, SCAD and MCP penalties, and
a coordinatewise root for any other separable parent with a smooth
log-density.

- [`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
  : Proximal Operator of a Penalty
- [`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
  : The Proximal Operator as a Piecewise Linear Table
- [`prox_apply()`](https://statmodels7.github.io/penalties7/reference/prox_apply.md)
  : Apply a Piecewise Linear Table
- [`has_prox()`](https://statmodels7.github.io/penalties7/reference/has_prox.md)
  : Does a Penalty Supply a Proximal Operator?

## Validation

Numerical checks of the derivatives of a penalty, meant above all for a
penalty written outside the package.

- [`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
  : Check a Penalty Numerically
