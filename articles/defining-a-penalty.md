# Defining a penalty

A penalty is $`\rho(D\beta;\theta)`$: a linear map, a scalar function of
the mapped coefficients, and hyperparameters of its own. In many R
packages that fit a penalized model, the penalty is written inside the
fitting function as a term of the objective, and its kinks, its rank and
the dependence of its value on the hyperparameters are not available
outside that function. In penalties7 the penalty is an object, and a
fitting layer reads these quantities from it to choose a fitting scheme.

This vignette builds a penalty that the package does not ship, writes
the methods that it has to supply, and shows what the base class and
[`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
provide in return. The example is the Berhu, the reverse Huber penalty
of Owen (2007):

``` math
\rho(\beta) = \begin{cases}
  \lambda\lvert\beta\rvert, & \lvert\beta\rvert \le \delta,\\[2pt]
  \lambda\dfrac{\beta^2 + \delta^2}{2\delta}, & \lvert\beta\rvert > \delta,
\end{cases}
```

which is the lasso near the origin and a ridge away from it. The kink at
zero sets small coefficients exactly to zero. In the tails the penalty
grows quadratically, so a large coefficient is shrunk in proportion to
its size, as under a ridge penalty, while the lasso shrinks every
nonzero coefficient by the same amount.

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

`params_bounds` holds open intervals, and `link_params` holds, for each
hyperparameter, the link that maps its interval onto the unconstrained
scale on which an outer search moves. `params_smooth` records the
hyperparameters in which the value is differentiable. It is stored with
the object, and
[`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
does not read it: the validator places its grid from
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md).

## The methods

Each of the generics below that takes `theta` aligns it before it
dispatches. The generic reorders the list to `penalty@params`, strips
stray names off the values, checks that each value is a single number
inside `params_bounds`, and, where it also takes `beta`, checks that
`beta` has length `penalty@n_coef`. A method therefore receives a list
in the penalty’s own order and does not align it again. One helper reads
the list and records the branch in which each coefficient falls, and the
six methods for the value and the derivatives share it:

``` r

.berhu <- function(penalty, beta, theta) {
  list(lam = theta$lambda, del = theta$delta,
       a = abs(beta), big = abs(beta) > theta$delta)
}
```

The next three methods supply the value, the gradient in $`\beta`$ and
the Hessian in $`\beta`$:

``` r

S7::method(penalty_value, BerhuPenalty) <- function(penalty, beta, theta, ...) {
  p <- .berhu(penalty, beta, theta)
  sum(ifelse(p$big, p$lam * (beta^2 + p$del^2) / (2 * p$del), p$lam * p$a))
}

S7::method(penalty_gradient, BerhuPenalty) <- function(penalty, beta, theta, ...) {
  p <- .berhu(penalty, beta, theta)
  ifelse(p$big, p$lam * beta / p$del, p$lam * sign(beta))
}

S7::method(penalty_hessian, BerhuPenalty) <- function(penalty, beta, theta, ...) {
  p <- .berhu(penalty, beta, theta)
  diag(ifelse(p$big, p$lam / p$del, 0), length(beta))
}
```

Three more methods supply the blocks in the hyperparameters: the
gradient, the Hessian and the mixed block
$`\partial^2\rho/\partial\beta\,\partial\theta`$, which an outer search
uses to estimate the hyperparameters. Each of the three methods includes
the `scale` argument in its formals, because an S7 method must include
the named arguments of its generic and these three generics have one.
The methods return the derivatives on the parameter scale, and with
`scale = "link"` the generic carries them to the link scale through the
links in `link_params`.

``` r

S7::method(penalty_grad_theta, BerhuPenalty) <- function(
    penalty, beta, theta, scale = c("parameter", "link"), ...) {
  p <- .berhu(penalty, beta, theta)
  list(lambda = sum(ifelse(p$big, (beta^2 + p$del^2) / (2 * p$del), p$a)),
       delta  = sum(ifelse(p$big,
                           p$lam * (p$del^2 - beta^2) / (2 * p$del^2), 0)))
}

S7::method(penalty_hess_theta, BerhuPenalty) <- function(
    penalty, beta, theta, scale = c("parameter", "link"), ...) {
  p <- .berhu(penalty, beta, theta)
  list(lambda_lambda = 0,
       delta_delta   = sum(ifelse(p$big, p$lam * beta^2 / p$del^3, 0)),
       lambda_delta  = sum(ifelse(p$big,
                                  (p$del^2 - beta^2) / (2 * p$del^2), 0)))
}

S7::method(penalty_cross, BerhuPenalty) <- function(
    penalty, beta, theta, scale = c("parameter", "link"), ...) {
  p <- .berhu(penalty, beta, theta)
  list(lambda = ifelse(p$big, beta / p$del, sign(beta)),
       delta  = ifelse(p$big, -p$lam * beta / p$del^2, 0))
}
```

The last two methods are
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
and
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md).
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
returns the positions of the points at which some derivative of the
value is discontinuous, in the units of $`D\beta`$, and not their
number. The Berhu has three such points. At the origin the first
derivative jumps from $`-\lambda`$ to $`\lambda`$, and at $`\pm\delta`$
the first derivative is continuous while the second jumps from $`0`$ to
$`\lambda/\delta`$. SCAD has five, for the same two reasons. The points
at $`\pm\delta`$ move with the hyperparameter, so the method reads
`delta` from `theta`.
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
returns `FALSE`, because the value defined above carries no normalizing
constant and is therefore not the negative log-density of a prior.

``` r

S7::method(penalty_kinks, BerhuPenalty) <- function(penalty, theta, ...) {
  c(0, -theta$delta, theta$delta)
}
S7::method(is_proper, BerhuPenalty) <- function(penalty, ...) FALSE

pen <- berhu_penalty(5)
b <- c(-2.5, -0.4, 0, 0.7, 3.1)
th <- list(lambda = 1.3, delta = 1)

penalty_kinks(pen, th)
#> [1]  0 -1  1
is_proper(pen)
#> [1] FALSE
penalty_value(pen, b, th)
#> [1] 13.039
penalty_gradient(pen, b, th)
#> [1] -3.25 -1.30  0.00  1.30  4.03
diag(penalty_hessian(pen, b, th))
#> [1] 1.3 0.0 0.0 0.0 1.3
```

### The shapes of the hyperparameter blocks

Every hyperparameter block is a named list, not a matrix or a data
frame.
[`penalty_grad_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
returns one number per hyperparameter, keyed by `penalty@params`.
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
returns one number per pair of hyperparameters, keyed by the two names
joined by an underscore: the diagonal pairs come first, and then each
remaining pair follows in the order of `penalty@params` (here
`"lambda_lambda"`, `"delta_delta"` and `"lambda_delta"`).
[`penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
returns one vector of length `n_coef` per hyperparameter. The generics
read these names when they carry the blocks to the link scale. A shipped
penalty returns the same shapes:

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
compares each closed form with a numerical reference and reports one row
per comparison:

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
#> 9        hess_theta[lambda_delta] vs numDeriv 8.851854e-12     OK
```

The rows for the gradient in $`\beta`$ and in each hyperparameter
differentiate the value with `numDeriv`. The rows for the Hessian in
$`\beta`$, for every entry of the Hessian in the hyperparameters
(including the off-diagonal entry `lambda_delta`) and for the mixed
block differentiate the analytic first derivatives instead, with the
Richardson extrapolation of `numDeriv`, so that the value is never
differenced twice.

### A numerical reference is invalid at a kink

If a coordinate of $`\beta`$ sits exactly on a kink, the check of the
Hessian fails on correct code:

``` r

check_penalty(pen, beta = b, theta = th, verbose = FALSE)[1:2, ]
#>                                 check    max_error status
#> 1                gradient vs numDeriv 3.083081e-11     OK
#> 2 hessian vs numDeriv on the gradient 1.000000e+00 FAILED
```

The vector `b` contains an exact zero. The gradient jumps from
$`-\lambda`$ to $`\lambda`$ there, so a central difference of the
gradient across zero measures the jump and not the formula. When `beta`
is left `NULL`,
[`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
draws it from a fixed seed and shifts the draw until every coordinate of
$`D\beta`$ lies more than `0.05` from every point that
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
returns, and the first run passed because of this shift. Moving the zero
is enough:

``` r

check_penalty(pen, beta = c(-2.5, -0.4, 0.05, 0.7, 3.1),
              theta = th, verbose = FALSE)[1:2, ]
#>                                 check    max_error status
#> 1                gradient vs numDeriv 3.192694e-10     OK
#> 2 hessian vs numDeriv on the gradient 9.635369e-12     OK
```

For this reason
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
returns the positions of the kinks and not their number:
[`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
needs the positions to keep its grid away from them.

## Methods supplied by the base class

``` r

c(is_quadratic   = is_quadratic(pen),
  beta_quadratic = beta_quadratic(pen, th),
  has_prox       = has_prox(pen))
#>   is_quadratic beta_quadratic       has_prox 
#>          FALSE          FALSE          FALSE
penalty_draw(pen, th)
#> [1] NA NA NA NA NA
try(penalty_rank(pen))
#> Error : Only a quadratic penalty has a rank; see is_quadratic().
```

The base class `penalty` has methods for the generics that a branch may
leave out.
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
and
[`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md)
return `FALSE`,
[`penalty_readable()`](https://statmodels7.github.io/penalties7/reference/penalty_readable.md)
and
[`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
return `NULL`, and
[`penalty_draw()`](https://statmodels7.github.io/penalties7/reference/penalty_draw.md)
returns `NA` for every coordinate of a penalty that is not quadratic.
[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md),
[`penalty_rank()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md),
[`penalty_null_basis()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
and
[`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
signal an error, because a rank, a null basis and a log
pseudo-determinant are properties of a quadratic form and the Berhu has
no matrix.
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
signals an error, and so do
[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md),
[`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.md)
and
[`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.md),
the derivatives that a Laplace approximation to a marginal likelihood
needs, and the three generics for the derivatives of the Hessian in
$`\beta`$
([`penalty_dhessian_beta()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian_beta.md),
[`penalty_d2hessian_beta()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.md)
and
[`penalty_dhessian_beta_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.md)).
The hyperparameters of the Berhu as written here therefore cannot be
estimated by a marginal criterion.

[`has_prox()`](https://statmodels7.github.io/penalties7/reference/has_prox.md)
returns `FALSE` because the Berhu has no method for
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md).
statmodels7 routes a block by whether
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
returns a point. A smooth penalty joins the Newton-type system of the
other coefficients, and a penalty with a kink is fitted by a proximal
gradient iteration, which calls
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md),
or by a coordinate descent, which reads
[`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md).
A penalty with a kink and without a method for
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
can be evaluated and differentiated away from its kinks, but it cannot
be fitted by either scheme.

The proximal operator of the Berhu has a closed form, piecewise linear
in $`v`$, so a method for
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
can be registered in the same way as the methods above.
[`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
returns the same operator as a piecewise-linear table, which the
compiled coordinate descent reads. The method of the base class returns
`NULL`, and statmodels7 then fits the block by the proximal gradient
iteration.

## The normalizing constant is kept

A separable penalty built from a distributions7 family is the negative
log-density of that family, and the shipped ones keep the constant. The
lasso is a Laplace at zero in its rate parametrization, so its value
carries $`-\log(\lambda/2)`$ per coefficient:

``` r

pen2 <- berhu_penalty(2)
small <- c(-0.3, 0.2)
c(berhu = penalty_value(pen2, small, th),
  lasso = penalty_value(lasso_penalty(n_coef = 2), small,
                        list(lambda = 1.3)),
  constant = -2 * log(1.3 / 2))
#>     berhu     lasso  constant 
#> 0.6500000 1.5115658 0.8615658
```

Without the constant the hyperparameters cannot be estimated: the
objective falls monotonically as $`\lambda \to 0`$, and the degrees of
freedom of a $`t`$ prior cannot be estimated at all. The Berhu above has
no constant, because its value is not a normalized density, and for that
reason its
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
method returns `FALSE`. Its hyperparameters are set by the caller or
chosen by a criterion on prediction error, not by a marginal likelihood.

## Summary

- Minimum to define a penalty: a subclass of `penalty` and methods for
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
  generics declare one. Each returns a named list on the parameter
  scale, and the generic carries it to the link scale.
- [`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
  returns positions, so that a numerical reference can be kept away from
  them.
- [`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
  compares the gradients with `numDeriv` on the value and the second
  derivatives with `numDeriv` on the analytic gradients, including every
  entry of the hyperparameter Hessian. When it is given no $`\beta`$ it
  draws one away from every declared kink.
- [`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md),
  [`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md),
  [`penalty_readable()`](https://statmodels7.github.io/penalties7/reference/penalty_readable.md),
  [`penalty_draw()`](https://statmodels7.github.io/penalties7/reference/penalty_draw.md),
  [`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md),
  [`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md),
  [`penalty_rank()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md),
  [`penalty_null_basis()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md),
  [`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md),
  [`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md),
  [`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
  and its siblings have methods on the base class. The first five return
  `FALSE`, `NULL` or `NA`; the others signal an error for a penalty that
  is not quadratic or has no such method.
- [`has_prox()`](https://statmodels7.github.io/penalties7/reference/has_prox.md)
  is `FALSE` until a method for
  [`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
  is registered, and a fit of a penalty with a kink needs that method.
