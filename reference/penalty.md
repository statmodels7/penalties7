# S7 Base Class for Penalties

This class is the abstract parent of every penalty in this package. A
penalty is a scalar function of the coefficients, written \\\rho(D\beta;
\theta)\\: a linear map \\D\\ that selects or combines coefficients, a
scalar function \\\rho\\ applied to what the map returns, and
hyperparameters \\\theta\\ that scale or shape it. The class holds no
mathematics of its own. It records the pieces that every branch needs,
so that a penalty's size, its hyperparameter names and their bounds can
be read in the same way for every branch.

## Usage

``` r
penalty(
  penalty_name = character(0),
  map = NULL,
  n_coef = integer(0),
  params = character(0),
  params_bounds = list(),
  link_params = list(),
  params_smooth = logical(0)
)
```

## Arguments

- penalty_name:

  A single string naming the penalty, used by
  [`print()`](https://rdrr.io/r/base/print.html) and by any code that
  reports which penalty a block carries.

- map:

  The matrix \\D\\, of \\m\\ rows and \\q\\ columns, or `NULL` for the
  identity. A Matrix object is kept in its own storage, so a diagonal
  map, such as the map that standardization produces, is recognized as
  diagonal by its class.

- n_coef:

  The number of coefficients \\q\\. The value is a single whole number.

- params:

  The hyperparameter names, in the order by which every derivative list
  is keyed. The value is `character(0)` for a penalty with none.

- params_bounds:

  A named list with one entry per hyperparameter, and each entry is a
  numeric pair giving an open interval. A value at either endpoint is
  rejected, so `(0, Inf)` excludes zero.

- link_params:

  A named list, one linkfunctions7 link per hyperparameter, carrying
  that hyperparameter's own interval onto the whole real line.

- params_smooth:

  A logical vector, one entry per hyperparameter, `TRUE` where the value
  is differentiable in it.

## Value

An S7 object of class `penalty` carrying the seven properties above. The
class is abstract: an object of exactly this class can be printed, it
receives the defaults listed under Methods, and
[`penalty_value()`](https://statmodels7.github.io/penalties7/reference/penalty_value.md)
and the other generics without a base method signal an error for it.

## The properties

With \\q\\ coefficients and a map of \\m\\ rows, `beta` is a vector of
length \\q\\, \\D\beta\\ has length \\m\\, and every derivative in the
coefficients is returned with length \\q\\ or shape \\q \times q\\. A
`NULL` map is the identity, and then \\m = q\\ and applying the map
involves no arithmetic.

Each hyperparameter carries a linkfunctions7 link in `link_params`,
which maps its own open interval onto the whole real line. A penalty can
therefore be optimized on the unconstrained scale: a caller works in
\\\eta = g(\theta)\\, without bound constraints, and
[`penalty_grad_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md),
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
and
[`penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
return either scale through their `scale` argument.

`params_smooth` records in which hyperparameters the value is
differentiable. It is `TRUE` for every hyperparameter of the quadratic,
structured, additive, SCAD and MCP branches. A separable penalty copies
the entries of its parent distribution, so a penalty built by
[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md)
on a Laplace distribution with a free location has `FALSE` for that
location.

## Writing a branch

This class is constructed directly only to write a new branch. A
subclass supplies
[`penalty_value()`](https://statmodels7.github.io/penalties7/reference/penalty_value.md),
[`penalty_gradient()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md),
[`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md),
[`penalty_grad_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md),
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md),
[`penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md),
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md)
and
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md),
which have no method on the base class, and
[`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
then compares the derivatives with numerical derivatives of the value.
The methods registered on the base class, listed below, supply defaults
for everything else.

The shipped branches are
[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
(with
[`ridge_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md)),
[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md)
(with
[`lasso_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md)
and
[`elasticnet_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md)),
[`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md),
[`mcp_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md),
[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
and
[`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md),
which is built on a parameters7 matrix parameter.

## Methods

The methods registered on `penalty` itself are inherited by every branch
that does not register its own.

Predicates and readable quantities:

- [`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
  returns `FALSE`.

- [`penalty_readable()`](https://statmodels7.github.io/penalties7/reference/penalty_readable.md)
  returns `NULL`.

- [`beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.penalty.md)
  returns `FALSE`.

Quantities of a quadratic penalty, which signal an error:

- [`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md),
  [`penalty_rank()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md),
  [`penalty_null_basis()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
  and
  [`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md).

Derivatives for a marginal criterion, which signal an error:

- [`penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.penalty.md),
  [`penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.penalty.md)
  and
  [`penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.penalty.md);

- [`penalty_dhessian_beta()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian_beta.penalty.md);

- [`penalty_d2hessian_beta()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.penalty.md)
  and
  [`penalty_dhessian_beta_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.penalty.md).

Proximal operator:

- [`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
  signals an error.

- [`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
  returns `NULL`.

Draws and printing:

- [`penalty_draw()`](https://statmodels7.github.io/penalties7/reference/penalty_draw.md)
  draws from a quadratic penalty read as a Gaussian prior and returns
  `NA` otherwise.

- [`print()`](https://statmodels7.github.io/penalties7/reference/print.penalty.md)
  writes one line describing the penalty.

## See also

[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md),
[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md),
[`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md),
[`mcp_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md),
[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
and
[`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md)
for the branches;
[`penalty_value()`](https://statmodels7.github.io/penalties7/reference/penalty_value.md)
and
[`penalty_gradient()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
for what a branch supplies;
[`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
to verify one.

## Examples

``` r
# Every branch inherits from this class.
pen <- quadratic_penalty(crossprod(diff(diag(4))), map = NULL)
S7::S7_inherits(pen, penalty)
#> [1] TRUE

# These properties can be read in the same way for every branch.
pen@penalty_name
#> [1] "quadratic"
pen@n_coef
#> [1] 4
pen@params
#> [1] "lambda"
pen@params_bounds
#> $lambda
#> [1]   0 Inf
#> 

# A first-difference map over three coefficients takes three coefficients
# in and returns two rows.
D <- diff(diag(3))
mapped <- quadratic_penalty(diag(2), map = D)
mapped
#> quadratic penalty on 3 coefficient(s) through 2 row(s); theta: lambda
```
