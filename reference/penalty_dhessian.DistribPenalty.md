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

The second derivatives read the parent's
[`distrib_grad_y_hess`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.html)
and
[`distrib_hess_y_hess`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.html),
so nothing is differentiated here either. A parent with closed forms for
those – the gaussian, hence every ridge and every random effect – makes
this branch exact; one without them inherits that package's documented
fallback, which is one central difference of its analytic first-order
component.
