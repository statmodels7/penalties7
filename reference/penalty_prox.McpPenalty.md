# Proximal Operator of MCP

The closed piecewise operator: a rescaled soft threshold below
\\\gamma\lambda\\ and the identity beyond it. Convex only while \\t \<
\gamma\\, and a longer step is rejected.

## Arguments

- pen:

  An `McpPenalty` object.

- v:

  A numeric vector.

- step:

  The step length.

- theta:

  A named list containing `lambda` and `gamma`.

- ...:

  Unused.

## Value

A numeric vector.
