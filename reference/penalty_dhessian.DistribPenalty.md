# Marginal Derivatives of a Separable Penalty

With \\\rho = -\sum_j \log f((D\beta)\_j;\theta)\\ the Hessian is
\\-D'\mathrm{diag}(\ell^{(yy)})D\\, so its \\\theta\\-derivatives are
the parent's `distrib_cross2_y` carried through the same map, and the
mixed third derivative is `distrib_cross_y` differentiated once more in
\\\theta\\. Both come from distributions7 rather than being written
again here.

## Arguments

- pen:

  A `DistribPenalty` object.

- beta:

  A numeric vector of coefficients.

- theta:

  A named list of hyperparameters.

- ...:

  Unused.

## Value

A named list of matrices, of vectors, or a logical.

## Details

The second derivatives use one central difference of the parent's
analytic `distrib_cross2_y` and `distrib_cross_y` in each
hyperparameter, the second \\\theta\\-derivative of a response
derivative not being one of that package's generics. A family that gains
it will be exact here with no change: the difference is taken of
whatever the parent supplies.
