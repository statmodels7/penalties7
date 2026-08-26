# Coefficient Derivatives of a SCAD Penalty

[`penalty_gradient()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
returns SCAD's defining derivative carried back through the map, and
[`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
its second derivative, which is \\-1/(a-1)\\ in the tapering region and
zero elsewhere.

## Arguments

- pen:

  A
  [`ScadPenalty()`](https://statmodels7.github.io/penalties7/reference/ScadPenalty.md)
  object.

- beta:

  A numeric vector of length `pen@n_coef`.

- theta:

  A list holding `lambda` and `a`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`penalty_gradient()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
a numeric vector of length `pen@n_coef`;
[`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
a symmetric `pen@n_coef` by `pen@n_coef` base matrix, diagonal under the
identity map and negative semidefinite always.

## Details

On each coordinate, with \\u = \lvert t\rvert\\ and \\s\\ its sign,

\$\$\frac{\partial\rho}{\partial t} = s \cdot \begin{cases} \lambda & u
\le \lambda \\ (a\lambda - u)/(a-1) & \lambda \< u \le a\lambda \\ 0 & u
\> a\lambda \end{cases}, \qquad \frac{\partial^2\rho}{\partial t^2} =
\begin{cases} -1/(a-1) & \lambda \< u \le a\lambda \\ 0 &
\text{otherwise.}\end{cases}\$\$

The Hessian is **negative** where it is not zero, which is the
non-convexity the family is named for: the penalty bends the wrong way
over the taper. A penalized objective stays convex only while the
likelihood's own curvature exceeds \\\lambda/(a-1)\\, and a solver that
assumes convexity has no guarantee here.

At \\t = 0\\ the sign is zero, so the gradient returned is `0`: one
element of the subdifferential \\\[-\lambda, \lambda\]\\, and the one a
smooth method can take a step with. Selection happens through
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md),
not through this.

## See also

[`penalty_value.ScadPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.ScadPenalty.md)
for the quantity differentiated,
[`penalty_grad_theta.ScadPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.ScadPenalty.md)
for the hyperparameter blocks,
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
for the step that produces exact zeros.

## Examples

``` r
pen <- scad_penalty(n_coef = 3)
th <- list(lambda = 1, a = 3.7)

# Linear, tapering, flat: the gradient falls to zero as u grows.
penalty_gradient(pen, c(0.5, 2, 5), th)
#> [1] 1.0000000 0.6296296 0.0000000

# The Hessian is -1/(a - 1) on the tapering coordinate alone.
diag(penalty_hessian(pen, c(0.5, 2, 5), th))
#> [1]  0.0000000 -0.3703704  0.0000000
-1 / (3.7 - 1)
#> [1] -0.3703704
```
