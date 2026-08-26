# Carry Hyperparameter Derivatives Onto the Unconstrained Scale

Applies the chain rule that turns a derivative in \\\theta\\ into one in
\\\eta = g(\theta)\\, at first and second order, using the diagonal
Jacobian the links supply. Handles one of the three derivative kinds per
call, whichever of `g`, `H` and `cross` is given.

## Usage

``` r
ptheta_to_link(pen, theta, g = NULL, H = NULL, cross = NULL)
```

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- theta:

  The aligned hyperparameters, as
  [`align_ptheta()`](https://statmodels7.github.io/penalties7/reference/align_ptheta.md)
  returns them.

- g:

  The parameter-scale gradient list, keyed by `pen@params`. Required for
  the gradient and for the Hessian; `NULL` otherwise.

- H:

  The parameter-scale Hessian list, keyed as
  [`ptheta_pairs()`](https://statmodels7.github.io/penalties7/reference/ptheta_pairs.md)
  keys it. Supply it with `g` to get the second order.

- cross:

  The parameter-scale mixed list, keyed by `pen@params`, each element a
  vector of length `pen@n_coef`. Supplied alone.

## Value

Whichever kind was supplied, on the unconstrained scale, with the same
names and shapes as the input. `cross` takes precedence over `H`, and
`H` over `g`, when more than one is given.

## Details

The links are scalar and one per hyperparameter, so the Jacobian is
diagonal and the chain rule needs no partition sums. Writing \\h =
g^{-1}\\ and \\\eta_i = g_i(\theta_i)\\,

\$\$\frac{\partial\rho}{\partial\eta_i} =
\frac{\partial\rho}{\partial\theta_i}\\ h_i'(\eta_i), \qquad
\frac{\partial^2\rho}{\partial\eta_i \partial\eta_j} =
\frac{\partial^2\rho}{\partial\theta_i \partial\theta_j}\\
h_i'(\eta_i)\\ h_j'(\eta_j) + \delta\_{ij}\\
\frac{\partial\rho}{\partial\theta_i}\\ h_i''(\eta_i).\$\$

The second-derivative term appears on the diagonal alone, which is why
the Hessian branch needs the gradient as well. The mixed block
\\\partial^2\rho / \partial\beta\\\partial\theta_i\\ is first order in
\\\theta\\ and picks up one factor of \\h_i'\\, the coefficient
direction being untouched by a reparametrization of the hyperparameters.

This is distributions7's interception, restricted to the two orders a
penalty consumer needs.
