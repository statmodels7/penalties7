# Hyperparameter Derivatives of a Penalty

The three blocks a consumer needs to move the hyperparameters rather
than the coefficients. `penalty_grad_theta()` returns
\\\partial\rho/\partial\theta_k\\, `penalty_hess_theta()` returns
\\\partial^2\rho/\partial\theta_k\partial\theta_l\\, and
`penalty_cross()` returns the mixed block
\\\partial^2\rho/\partial\beta\\\partial\theta_k\\, one coefficient
vector per hyperparameter. Each answers on the parameter scale or on the
unconstrained scale, chosen by `scale`.

## Usage

``` r
penalty_grad_theta(pen, beta, theta, scale = c("parameter", "link"), ...)

penalty_hess_theta(pen, beta, theta, scale = c("parameter", "link"), ...)

penalty_cross(pen, beta, theta, scale = c("parameter", "link"), ...)
```

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object of any branch.

- beta:

  A numeric vector of length `pen@n_coef`, coerced with
  [`as.numeric()`](https://rdrr.io/r/base/numeric.html) in the generic.

- theta:

  A named list of hyperparameter values, or a named numeric vector
  carrying the same, holding every name in `pen@params`.

- scale:

  `"parameter"`, the default, gives derivatives in \\\theta\\. `"link"`
  gives them in the unconstrained values \\\eta = g(\theta)\\, which is
  the scale an optimizer works on and the scale a hyperparameter's
  standard error is built on.

- ...:

  Passed to methods. No shipped method reads it.

## Value

All three return a named list. `penalty_grad_theta()` has one element
per hyperparameter, each a single number, named by `pen@params`.
`penalty_hess_theta()` has \\p(p+1)/2\\ elements, each a single number,
keyed by pairs as described above. `penalty_cross()` has one element per
hyperparameter, each a numeric vector of length `pen@n_coef`.

## What each is for

The gradient and the Hessian in \\\theta\\ are what a joint maximization
of coefficients and hyperparameters steps with, and what the Laplace
term of a marginal criterion differentiates. The mixed block is the
off-diagonal corner of that joint system, and it is also what the
implicit function theorem needs: the penalized mode moves with a
hyperparameter as \\\partial\hat\beta/\partial\theta_k = -(H + S)^{-1}
\partial^2\rho/\partial\beta\\\partial\theta_k\\, so a criterion
computed at the mode cannot be differentiated without it.

## The two scales

Each hyperparameter carries a link \\g\\ mapping its open interval onto
the whole line. With `scale = "link"` the derivatives are taken in
\\\eta = g(\theta)\\, the coordinates an unconstrained optimizer moves.
The chain rule is applied in the generic body, so a method always
returns the parameter scale and a branch written by a user gets both
scales for free. Writing \\h = g^{-1}\\,

\$\$\frac{\partial\rho}{\partial\eta_k} =
\frac{\partial\rho}{\partial\theta_k} h_k'(\eta_k), \qquad
\frac{\partial^2\rho}{\partial\eta_k \partial\eta_l} =
\frac{\partial^2\rho}{\partial\theta_k \partial\theta_l} h_k'(\eta_k)
h_l'(\eta_l) + \delta\_{kl} \frac{\partial\rho}{\partial\theta_k}
h_k''(\eta_k),\$\$

with the mixed block picking up one factor of \\h_k'\\ and no
second-order term, a reparametrization of the hyperparameters leaving
the coefficient direction alone.

## The keying

The gradient and the mixed block are keyed by `pen@params`. The Hessian
is keyed by pairs, **diagonals first** in `params` order and then the
upper off-diagonal pairs, each joined by an underscore: for an elastic
net the keys are `lambda_lambda`, `alpha_alpha`, `lambda_alpha`. Only
the upper triangle is returned, the matrix being symmetric.

## See also

[`penalty_value()`](https://statmodels7.github.io/penalties7/reference/penalty_value.md)
for the quantity differentiated,
[`penalty_gradient()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
for the coefficient derivatives,
[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
for the third-order blocks a marginal criterion reads,
[`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
to verify all three against numerical routes.

## Examples

``` r
pen <- quadratic_penalty(diag(2))
penalty_grad_theta(pen, c(1, -1), list(lambda = 3))
#> $lambda
#> [1] 0.6666667
#> 
penalty_hess_theta(pen, c(1, -1), list(lambda = 3))
#> $lambda_lambda
#> [1] 0.1111111
#> 
penalty_cross(pen, c(1, -1), list(lambda = 3))
#> $lambda
#> [1]  1 -1
#> 

# An elastic net has two hyperparameters, so three Hessian keys.
enet <- elasticnet_penalty(n_coef = 3)
names(penalty_hess_theta(enet, c(1, -0.4, 0.2),
                         list(lambda = 2, alpha = 0.4)))
#> [1] "lambda_lambda" "alpha_alpha"   "lambda_alpha" 

# The link scale is the parameter scale times the chain factor. lambda
# rides a log link, so h'(eta) is lambda itself.
g_par  <- penalty_grad_theta(pen, c(1, -1), list(lambda = 3))
g_link <- penalty_grad_theta(pen, c(1, -1), list(lambda = 3),
                             scale = "link")
all.equal(g_link$lambda, g_par$lambda * 3)
#> [1] TRUE
```
