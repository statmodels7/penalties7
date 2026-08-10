# Marginal Derivatives of an Additive Penalty

\\S = \sum_k \lambda_k P_k\\ is linear in the smoothing parameters, so
\\\partial S/\partial\lambda_k = P_k\\ and every higher derivative is
zero.

## Arguments

- pen:

  An `AdditivePenalty` object.

- beta:

  A numeric vector of coefficients.

- theta:

  A named list of hyperparameters.

- ...:

  Unused.

## Value

A named list of matrices, of vectors, or a logical.
