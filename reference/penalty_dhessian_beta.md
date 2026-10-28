# The Derivative of the Coefficient Hessian in the Coefficients

Returns \\\partial S/\partial\beta\\\[v\] = \sum_c v_c\\
\partial^3\rho/\partial\beta\\\partial\beta'\\\partial\beta_c\\, the
derivative of the coefficient Hessian along a direction \\v\\, as one
matrix of side `penalty@n_coef`. It is zero for every penalty whose
Hessian does not depend on the coefficients, apart from a kinked
separable penalty, for which the generic signals an error.

## Usage

``` r
penalty_dhessian_beta(penalty, beta, theta, v, ...)
```

## Arguments

- penalty:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- beta:

  A numeric vector of length `penalty@n_coef`; the generic signals an
  error for any other length.

- theta:

  A named list of hyperparameter values, or a named numeric vector
  carrying the same, checked by the generic as in
  [`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md).

- v:

  A numeric vector of length `penalty@n_coef`, the direction; the
  generic signals an error for any other length.

- ...:

  Passed to methods. The shipped methods do not read it.

## Value

A square base matrix of side `penalty@n_coef`.

## Role in a marginal criterion

A marginal criterion reads \\\log\lvert H + S\rvert\\ at the penalized
mode, and the mode moves with the hyperparameters. When \\S\\ does not
depend on the coefficients, that movement reaches the determinant
through \\H\\ alone. When it does, as for a heavy-tailed prior on a
random effect, the movement reaches the determinant through \\S\\ as
well, and the gradient contains \\\mathrm{tr}(M\\\partial
S/\partial\beta\[v\])\\ with \\v\\ the derivative of the mode. A
prediction-error criterion contains the same matrix inside the trace
that gives its degrees of freedom.

The method contracts the direction and does not return the array of
order three, which has \\q^3\\ entries, because the criteria read that
array only along a direction.

## The value on each branch

|  |  |
|----|----|
| branch | \\\partial S/\partial\beta\[v\]\\ |
| [`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md), [`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md), [`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md) | zero |
| [`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md) whose parent is quadratic in the argument | zero |
| [`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md) with a univariate parent otherwise | \\-D'\mathrm{diag}(\ell^{(yyy)}(D\beta)\odot Dv)\\D\\ |
| [`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md) with a multivariate parent otherwise | \\-D'\mathrm{blockdiag}(T_i)D\\, \\(T_i)\_{ab} = \sum_c \ell^{(y_ay_by_c)}(b_i)(Dv)\_{ic}\\ |
| a kinked parent, [`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md), [`mcp_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md) | an error |

The univariate row follows from \\S =
-D'\mathrm{diag}(\ell^{(yy)}(D\beta))D\\: differentiating
\\\ell^{(yy)}((D\beta)\_j)\\ in \\\beta_c\\ gives
\\\ell^{(yyy)}((D\beta)\_j) D\_{jc}\\, and summing against \\v_c\\ gives
the factor \\(Dv)\_j\\. The third response derivative is read from
[`distributions7::distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.html),
which is closed form for every family that distributions7 ships, the
Student t among them, and one stencil on the log-density for a family
defined by its density alone.

In the multivariate row \\b_i\\ is block \\i\\ of \\D\beta\\ and
\\(Dv)\_{ic}\\ is coordinate \\c\\ of block \\i\\ of \\Dv\\. A
multivariate parent that is not quadratic, such as a multivariate
Student t prior, is read through
[`distributions7::distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.html),
which returns one \\p \times p \times p\\ array per block;
[`dp_contract()`](https://statmodels7.github.io/penalties7/reference/dp_contract.md)
contracts its last index against the coordinates of the block in \\Dv\\.
If the parent cannot supply the array, the method signals an error that
names the penalty (see
[`dp_parent_or_reject()`](https://statmodels7.github.io/penalties7/reference/dp_contract.md)).
The separable method calls
[`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md)
first, which evaluates the response Hessian of the parent, and returns
zero when it gives `TRUE`. A Gaussian prior of any dimension therefore
gives zero, and the third response derivative of the parent is not
computed.

## See also

[`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
for the quantity differentiated,
[`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.md)
for the predicate under which it is zero,
[`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.md)
for the derivative in the hyperparameters,
[`penalty_d2hessian_beta()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.md)
for the derivatives of this quantity.

## Examples

``` r
b <- c(1, -0.5, 0.3)
v <- c(0.2, 0.1, -0.4)

# The derivative is zero for a quadratic penalty, whose Hessian does not
# change with beta.
penalty_dhessian_beta(quadratic_penalty(diag(3)), b, list(lambda = 2), v)
#>      [,1] [,2] [,3]
#> [1,]    0    0    0
#> [2,]    0    0    0
#> [3,]    0    0    0

# The Hessian of a Student t prior changes with beta, and this is its
# derivative along v, compared with a central difference.
h <- heavy_penalty(n_coef = 3)
th <- list(sigma = 1, nu = 4)
D <- penalty_dhessian_beta(h, b, th, v)
eps <- 1e-6
num <- (penalty_hessian(h, b + eps * v, th) -
        penalty_hessian(h, b - eps * v, th)) / (2 * eps)
max(abs(D - num))
#> [1] 8.447351e-11

# A kinked parent has no such derivative at the kink, and the call is an
# error.
try(penalty_dhessian_beta(lasso_penalty(n_coef = 3), b, list(lambda = 1), v))
#> Error : 'separable [fixed laplace2 [mu=0]]' has a kink, so penalty_dhessian_beta() does not exist there and
#>   its hyperparameters cannot be estimated by a marginal criterion.
```
