# Defining a penalty

A penalty is $`\rho(D\beta;\theta)`$: a linear map, a scalar function of
the mapped coefficients, and hyperparameters of its own. Almost every R
package that fits a penalized model writes one inside the function that
fits it, so the penalty is a term in an objective rather than an object,
and nothing can ask it whether it has a kink, what its rank is, or how
its value moves with its hyperparameter. Those are the questions a
fitting layer has to answer before it can choose a scheme.

This vignette builds a penalty the package does not ship and shows what
it has to supply and what it gets in return. The example is the Berhu,
the reverse Huber of Owen (2007):

``` math
\rho(\beta) = \begin{cases}
  \lambda\lvert\beta\rvert, & \lvert\beta\rvert \le \delta,\\[2pt]
  \lambda\dfrac{\beta^2 + \delta^2}{2\delta}, & \lvert\beta\rvert > \delta,
\end{cases}
```

the lasso near the origin and the ridge away from it. It has a kink at
zero, so it selects; it grows quadratically in the tails, so a large
coefficient is not shrunk proportionally.

## The class and its constructor

``` r

BerhuPenalty <- S7::new_class("BerhuPenalty", parent = penalty)

berhu_penalty <- function(n_coef, map = NULL) {
  BerhuPenalty(
    penalty_name = "berhu",
    map = map,
    n_coef = as.integer(n_coef),
    params = c("lambda", "delta"),
    params_bounds = list(lambda = c(0, Inf), delta = c(0, Inf)),
    link_params = list(lambda = linkfunctions7::log_link(),
                       delta = linkfunctions7::log_link()),
    params_smooth = c(lambda = TRUE, delta = TRUE)
  )
}
```

`params_bounds` are open intervals and `link_params` says how each
hyperparameter reaches the unconstrained scale an outer search steps on.
`params_smooth` records which of them the value is differentiable in, so
a numerical guard can switch itself off where it would be measuring a
kink.

## The methods

Every generic aligns `theta` before it dispatches, reordering it to
`pen@params`, stripping stray names off the values and checking it
against `params_bounds`. A method therefore receives a list already in
the family’s own order and does not align it again. One helper reads it
and works out the branch each coefficient falls in, so the eight methods
share it:

``` r

.berhu <- function(pen, beta, theta) {
  list(lam = theta$lambda, del = theta$delta,
       a = abs(beta), big = abs(beta) > theta$delta)
}
```

The value, the gradient in $`\beta`$ and the Hessian in $`\beta`$:

``` r

S7::method(penalty_value, BerhuPenalty) <- function(pen, beta, theta, ...) {
  p <- .berhu(pen, beta, theta)
  sum(ifelse(p$big, p$lam * (beta^2 + p$del^2) / (2 * p$del), p$lam * p$a))
}

S7::method(penalty_gradient, BerhuPenalty) <- function(pen, beta, theta, ...) {
  p <- .berhu(pen, beta, theta)
  ifelse(p$big, p$lam * beta / p$del, p$lam * sign(beta))
}

S7::method(penalty_hessian, BerhuPenalty) <- function(pen, beta, theta, ...) {
  p <- .berhu(pen, beta, theta)
  diag(ifelse(p$big, p$lam / p$del, 0), length(beta))
}
```

Then the three blocks in the hyperparameters, which are what a marginal
criterion needs in order to estimate them. Note the `scale` argument: an
S7 method’s formals have to include the generic’s named arguments, and
these three generics carry one.

``` r

S7::method(penalty_grad_theta, BerhuPenalty) <- function(
    pen, beta, theta, scale = c("parameter", "link"), ...) {
  p <- .berhu(pen, beta, theta)
  list(lambda = sum(ifelse(p$big, (beta^2 + p$del^2) / (2 * p$del), p$a)),
       delta  = sum(ifelse(p$big,
                           p$lam * (p$del^2 - beta^2) / (2 * p$del^2), 0)))
}

S7::method(penalty_hess_theta, BerhuPenalty) <- function(
    pen, beta, theta, scale = c("parameter", "link"), ...) {
  p <- .berhu(pen, beta, theta)
  list(lambda_lambda = 0,
       delta_delta   = sum(ifelse(p$big, p$lam * beta^2 / p$del^3, 0)),
       lambda_delta  = sum(ifelse(p$big,
                                  (p$del^2 - beta^2) / (2 * p$del^2), 0)))
}

S7::method(penalty_cross, BerhuPenalty) <- function(
    pen, beta, theta, scale = c("parameter", "link"), ...) {
  p <- .berhu(pen, beta, theta)
  list(lambda = ifelse(p$big, beta / p$del, sign(beta)),
       delta  = ifelse(p$big, -p$lam * beta / p$del^2, 0))
}
```

Finally the two structural questions.
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
returns the **positions** of the kinks, not their number; the Berhu’s
only kink is the origin, where SCAD has five.

``` r

S7::method(penalty_kinks, BerhuPenalty) <- function(pen, theta, ...) 0
S7::method(is_proper, BerhuPenalty) <- function(pen) TRUE

pen <- berhu_penalty(5)
b <- c(-2.5, -0.4, 0, 0.7, 3.1)
th <- list(lambda = 1.3, delta = 1)

penalty_value(pen, b, th)
#> [1] 13.039
penalty_gradient(pen, b, th)
#> [1] -3.25 -1.30  0.00  1.30  4.03
diag(penalty_hessian(pen, b, th))
#> [1] 1.3 0.0 0.0 0.0 1.3
```

### The shapes are part of the contract

Every hyperparameter block is a **named list**, not a matrix or a data
frame:
[`penalty_grad_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
one scalar per hyperparameter,
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
keyed `"lambda_lambda"`, `"delta_delta"`, `"lambda_delta"`, and
[`penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
one vector of length `n_coef` per hyperparameter. Compare a shipped
penalty and the shapes are the same:

``` r

str(penalty_hess_theta(scad_penalty(n_coef = 3), c(-2, 0.3, 1.5),
                       list(lambda = 1, a = 3.7)))
#> List of 3
#>  $ lambda_lambda: num -0.741
#>  $ a_a          : num -0.0635
#>  $ lambda_a     : num -0.206
```

## Validating it

[`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
differentiates the value numerically and compares, one block at a time:

``` r

print(check_penalty(pen, theta = th, verbose = FALSE))
#>                                         check    max_error status
#> 1                        gradient vs numDeriv 1.443589e-11     OK
#> 2         hessian vs numDeriv on the gradient 8.419077e-12     OK
#> 3              grad_theta[lambda] vs numDeriv 5.702300e-12     OK
#> 4       hess_theta[lambda_lambda] vs numDeriv 0.000000e+00     OK
#> 5 cross[lambda] vs Richardson on the gradient 2.740665e-12     OK
#> 6               grad_theta[delta] vs numDeriv 4.268279e-13     OK
#> 7         hess_theta[delta_delta] vs numDeriv 9.007713e-12     OK
#> 8  cross[delta] vs Richardson on the gradient 7.601645e-12     OK
```

The reference is `numDeriv`, which shares no arithmetic with the
formulas, and the cross block is checked by Richardson extrapolation on
the analytic gradient, never by differencing the value twice.

### A numerical reference is invalid at a kink

Hand the validator a $`\beta`$ with a coordinate sitting exactly on the
kink and the second-derivative check fails on correct code:

``` r

check_penalty(pen, beta = b, theta = th, verbose = FALSE)[1:2, ]
#>                                 check    max_error status
#> 1                gradient vs numDeriv 3.083081e-11     OK
#> 2 hessian vs numDeriv on the gradient 1.000000e+00 FAILED
```

`b` contains an exact zero, the gradient jumps from $`-\lambda`$ to
$`+\lambda`$ there, and a central difference across the jump measures
the kink, not the formula. When `beta` is left `NULL` the validator
draws it and nudges each draw off any position
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
declares, which is why the run above passed. Moving the zero is enough:

``` r

check_penalty(pen, beta = c(-2.5, -0.4, 0.05, 0.7, 3.1),
              theta = th, verbose = FALSE)[1:2, ]
#>                                 check    max_error status
#> 1                gradient vs numDeriv 3.192694e-10     OK
#> 2 hessian vs numDeriv on the gradient 9.635369e-12     OK
```

This is the reason
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
returns positions. A count would say a kink exists; the positions say
where a reference must not be taken.

## What the base class supplies, and what it refuses

``` r

c(is_quadratic   = is_quadratic(pen),
  beta_quadratic = beta_quadratic(pen),
  has_prox       = has_prox(pen))
#>   is_quadratic beta_quadratic       has_prox 
#>          FALSE          FALSE          FALSE
try(penalty_rank(pen))
#> Error : Only a quadratic penalty has a rank; see is_quadratic().
```

[`penalty_rank()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
and
[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
are refused because there is no matrix: a rank is a property of a
quadratic form, and returning one for a penalty that is not quadratic
would be a different object under the same name.

[`has_prox()`](https://statmodels7.github.io/penalties7/reference/has_prox.md)
is `FALSE` because no proximal operator has been registered, and that is
the answer with consequences. A fitting layer routes a block by whether
its penalty is twice differentiable: a smooth one joins the system
solved by `iwls()` or an ordinary optimizer, and a kinked one needs a
proximal step or a coordinate descent, both of which need the operator.
A kinked penalty with no
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
can be evaluated and differentiated away from its kinks, and cannot be
fitted by the scheme its kink calls for.

The Berhu’s operator is closed, so registering one is a further
exercise, not an obstruction;
[`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
is the piecewise-linear table a compiled coordinate descent reads, and
supplying it is what puts a penalty on that path.

## The normalizing constant is kept

A separable penalty built from a distributions7 family is the negative
log-density of that family, and the shipped ones keep the constant. The
lasso is a fixed Laplace in its rate chart, so its value carries
$`-\log(\lambda/2)`$ per coefficient:

``` r

small <- c(-0.3, 0.2)
c(berhu = penalty_value(pen, small, th),
  lasso = penalty_value(lasso_penalty(n_coef = 2), small,
                        list(lambda = 1.3)),
  constant = -2 * log(1.3 / 2))
#>     berhu     lasso  constant 
#> 0.6500000 1.5115658 0.8615658
```

Dropping the constant makes hyperparameter estimation degenerate:
without it the objective falls monotonically as $`\lambda \to 0`$, and a
$`t`$ prior’s degrees of freedom cannot be estimated at all. The Berhu
above has no constant because it is not a normalized density; a penalty
of that kind states its hyperparameters as what they are, quantities the
caller sets or an outer criterion chooses on a prediction error instead
of a marginal likelihood.

## Summary

- **Minimum to define a penalty:** a subclass of `penalty` and methods
  for
  [`penalty_value()`](https://statmodels7.github.io/penalties7/reference/penalty_value.md),
  [`penalty_gradient()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md),
  [`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md),
  [`penalty_grad_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md),
  [`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md),
  [`penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md),
  [`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
  and
  [`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md).
- The three hyperparameter blocks take a `scale` argument, because their
  generics declare one, and every one returns a **named list**.
- [`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
  returns positions, so a numerical reference can be kept off them; a
  count would not be enough.
- [`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
  compares against `numDeriv` and Richardson, and draws a $`\beta`$
  nudged off every declared kink when it is given none.
- [`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md),
  [`penalty_rank()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md),
  [`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md),
  [`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md),
  [`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
  and
  [`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
  come from the base class, and the first three refuse where the penalty
  is not quadratic.
- [`has_prox()`](https://statmodels7.github.io/penalties7/reference/has_prox.md)
  decides which fitting scheme a kinked block can reach.
