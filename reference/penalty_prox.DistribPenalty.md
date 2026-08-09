# Proximal Operator of a Separable Penalty

Closed form for the Gaussian (a shrinkage), the Laplace (a soft
threshold) and the elastic net (the one followed by the other); for any
other parent with a differentiable log-density, the coordinatewise root
of \\(\beta - v)/t = \ell^{(y)}(\beta)\\, which is unique because the
density is log-concave.

## Arguments

- pen:

  A `DistribPenalty` object.

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
