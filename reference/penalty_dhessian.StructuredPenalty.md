# Marginal Derivatives of a Structured Penalty

\\S = \Omega(\theta)\\ is the matrix parameter itself, so its
derivatives in the hyperparameters are the structure's own `param_d1`
and `param_d2`, which parameters7 supplies exactly. It is quadratic in
the coefficients, so the mixed third derivative is
\\\partial^2\Omega/\partial\theta_m\partial\theta_l\\\beta\\.

## Arguments

- pen:

  A `StructuredPenalty` object.

- beta:

  A numeric vector of coefficients.

- theta:

  A named list of hyperparameters.

- ...:

  Unused.

## Value

A named list of matrices, of vectors, or a logical.
