# Proximal Operator of a Quadratic Penalty

One linear solve, \\(I + tS)^{-1}v\\ with \\S = \lambda D'PD\\ the
penalty's Hessian, which holds for any map because the objective stays
quadratic.

## Arguments

- pen:

  A `QuadraticPenalty` object.

- v:

  A numeric vector.

- step:

  The step length.

- theta:

  A named list of hyperparameter values.

- ...:

  Unused.

## Value

A numeric vector.
