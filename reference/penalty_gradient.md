# Coefficient Derivatives of a Penalty

`penalty_gradient()` returns \\\partial\rho/\partial\beta\\ and
`penalty_hessian()` returns \\\partial^2\rho/\partial\beta^2\\, both in
closed form for every shipped branch and both already carried back
through the map. These are the two quantities a penalized fit adds to
the score and to the information at each iteration.

## Usage

``` r
penalty_gradient(pen, beta, theta, ...)

penalty_hessian(pen, beta, theta, ...)
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

- ...:

  Passed to methods. No shipped method reads it.

## Value

`penalty_gradient()` a numeric vector of length `pen@n_coef`;
`penalty_hessian()` a symmetric matrix of that side. It is a base matrix
even where the map is a Matrix, and a `dgCMatrix` for a penalty built by
[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
with `blocks > 1`, which keeps its sparse storage deliberately.

## Through the map

\\\rho\\ is a function of \\t = D\beta\\, so the chain rule gives

\$\$\frac{\partial\rho}{\partial\beta} = D^\top
\frac{\partial\rho}{\partial t}, \qquad
\frac{\partial^2\rho}{\partial\beta^2} = D^\top
\frac{\partial^2\rho}{\partial t^2} D,\$\$

and the map is applied here, once, so a caller adds the result straight
into a system in \\\beta\\. For a separable penalty the middle matrix is
diagonal and the product is formed without building it.

## At a kink

Where \\\rho\\ is not differentiable the returned gradient is one
element of the subdifferential, chosen by the branch, and the Hessian is
whatever the branch's own formula gives there. A lasso at \\\beta_j =
0\\ is the case that matters: the subdifferential is the whole interval
\\\[-\lambda, \lambda\]\\ and a single number cannot stand for it. Route
the block to
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
instead of to a gradient method, and use
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
to find out where the difficulty is.

## See also

[`penalty_value()`](https://statmodels7.github.io/penalties7/reference/penalty_value.md)
for the quantity differentiated,
[`penalty_grad_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
for the hyperparameter derivatives,
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
for the operator to use where the gradient does not exist,
[`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
to verify both against a numerical route.

## Examples

``` r
pen <- quadratic_penalty(diag(2))
penalty_gradient(pen, c(1, -1), list(lambda = 3))
#> [1]  3 -3
penalty_hessian(pen, c(1, -1), list(lambda = 3))
#>      [,1] [,2]
#> [1,]    3    0
#> [2,]    0    3

# A quadratic penalty has a constant Hessian, so the gradient is linear.
all.equal(penalty_gradient(pen, c(1, -1), list(lambda = 3)),
          drop(penalty_hessian(pen, c(1, -1), list(lambda = 3)) %*% c(1, -1)))
#> [1] TRUE

# The lasso's gradient is lambda times a sign, and at zero it is one
# element of the interval [-lambda, lambda] rather than a derivative.
penalty_gradient(lasso_penalty(n_coef = 3), c(2, -0.5, 0),
                 list(lambda = 1.5))
#> [1]  1.5 -1.5  0.0
```
