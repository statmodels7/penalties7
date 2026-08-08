# Proximal Operator of SCAD

The closed piecewise operator: a soft threshold near zero, a rescaled
threshold on the tapering region, and the identity beyond \\a\lambda\\.
The middle region is convex only while \\t \< a - 1\\, and a longer step
is rejected.

## Arguments

- pen:

  A `ScadPenalty` object.

- v:

  A numeric vector.

- step:

  The step length.

- theta:

  A named list containing `lambda` and `a`.

- ...:

  Unused.

## Value

A numeric vector.
