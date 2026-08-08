# Proximal Operator of a Penalty

The point that minimizes the penalty plus a quadratic pull towards
\\v\\, \$\$\mathrm{prox}\_{t\rho}(v) = \arg\min\_{\beta} \left\\
\tfrac{1}{2t}\lVert \beta - v \rVert^{2} + \rho(\beta;\theta)
\right\\,\$\$ which is what a proximal-gradient method applies after
each gradient step, and the only operation that lets a
non-differentiable penalty be minimized without differencing it.

## Usage

``` r
penalty_prox(pen, v, step, theta, ...)
```

## Arguments

- pen:

  A
  [`penalty`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- v:

  A numeric vector, the point to be pulled towards.

- step:

  The positive step length \\t\\.

- theta:

  A named list of hyperparameter values.

- ...:

  Passed to methods.

## Value

A numeric vector of the same length as `v`.

## Details

Where the operator has a closed form, that form is used. A quadratic
penalty gives one linear solve, \\(I + tS)^{-1}v\\ with \\S\\ the
penalty's Hessian, and works for any linear map. The separable branches
need the identity map, because with a general \\D\\ the problem does not
split by coordinate and has no elementary solution: that case is
rejected rather than approximated.

The named separable instances are exact: the Gaussian gives
\\v/(1+t/\sigma^{2})\\ and the Laplace the soft threshold
\\\mathrm{sign}(v)(\lvert v \rvert - t\lambda)\_{+}\\. SCAD and MCP have
the closed piecewise operators of their defining papers, each valid only
while the quadratic pull dominates the concave region of the penalty –
\\t \< a-1\\ for SCAD and \\t \< \gamma\\ for MCP – so a step beyond
that is rejected, the operator being set-valued there.

Any other separable penalty built from a distribution with a
differentiable log-density is solved coordinatewise from the stationary
condition \\(\beta - v)/t = \ell^{(y)}(\beta)\\. The right-hand side is
the response derivative
[`distrib_grad_y`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.html),
closed form for every continuous family, and log-concavity of the
density makes the left side minus the right side strictly increasing, so
the root is unique and a bracketed search finds it.

## See also

[`has_prox`](https://statmodels7.github.io/penalties7/reference/has_prox.md),
[`penalty_value`](https://statmodels7.github.io/penalties7/reference/penalty_value.md)

## Examples

``` r
# the lasso: a soft threshold
penalty_prox(lasso_penalty(n_coef = 3), c(2, 0.3, -1.4), 1, list(lambda = 1))
#> [1]  1.0  0.0 -0.4

# the ridge: a shrinkage
penalty_prox(ridge_penalty(n_coef = 2), c(2, -1), 1, list(sigma = 1))
#> [1]  1.0 -0.5
```
