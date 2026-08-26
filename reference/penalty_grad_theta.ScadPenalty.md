# Hyperparameter Derivatives of a SCAD Penalty

The three blocks in \\(\lambda, a)\\, each a closed piecewise form
summed over the coordinates.
[`penalty_grad_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
returns the two first derivatives,
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
the three second derivatives, and
[`penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
the mixed block, one coefficient vector per hyperparameter.

## Arguments

- pen:

  A
  [`ScadPenalty()`](https://statmodels7.github.io/penalties7/reference/ScadPenalty.md)
  object.

- beta:

  A numeric vector of length `pen@n_coef`.

- theta:

  A list holding `lambda` and `a`.

- scale:

  Read by the generic, which applies the chain rule onto the
  unconstrained scale after dispatch. These methods always return the
  parameter scale.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`penalty_grad_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
a list of two numbers, `lambda` and `a`.
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
a list of three numbers, `lambda_lambda`, `a_a` and `lambda_a`.
[`penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
a list of two numeric vectors of length `pen@n_coef`, named `lambda` and
`a`.

## Details

Each region contributes its own expression and the sums run over the
coordinates that fall in it. In \\u \le \lambda\\ the value is \\\lambda
u\\, so \\\partial\rho/\partial\lambda = u\\ and the shape does not
enter; beyond \\a\lambda\\ the value is \\(a+1)\lambda^2/2\\, so
\\\partial\rho/\partial\lambda = (a+1)\lambda\\ and
\\\partial\rho/\partial a = \lambda^2/2\\, neither depending on the
coefficient at all.

These derivatives exist, and no marginal criterion can use them: a
Laplace expansion needs the penalty to be twice differentiable at the
mode, and a mode with a coefficient held at zero sits on the kink. What
they serve is a joint step in which the hyperparameters are estimated
alongside the coefficients, and
[`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md),
which compares them against numDeriv.

## See also

[`penalty_value.ScadPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.ScadPenalty.md)
for the quantity differentiated,
[`penalty_grad_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
for the generic and its two scales,
[`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md),
which verifies all three.

## Examples

``` r
pen <- scad_penalty(n_coef = 3)
th <- list(lambda = 1, a = 3.7)
b <- c(0.5, 2, 5)

penalty_grad_theta(pen, b, th)
#> $lambda
#> [1] 7.57037
#> 
#> $a
#> [1] 0.5685871
#> 
penalty_hess_theta(pen, b, th)
#> $lambda_lambda
#> [1] 4.32963
#> 
#> $a_a
#> [1] -0.05080526
#> 
#> $lambda_a
#> [1] 0.8628258
#> 

# A coordinate past a * lambda contributes lambda^2 / 2 to the shape
# derivative and nothing to the mixed block.
penalty_cross(pen, b, th)$a
#> [1] 0.0000000 0.1371742 0.0000000
```
