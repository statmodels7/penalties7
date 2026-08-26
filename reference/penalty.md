# S7 Base Class for Penalties

The abstract parent of every penalty in this package. A penalty is a
scalar function of the coefficients, written \\\rho(D\beta; \theta)\\: a
linear map \\D\\ that selects or combines coefficients, a scalar
function \\\rho\\ applied to what the map returns, and hyperparameters
\\\theta\\ that scale or shape it. The class holds no mathematics of its
own; it records the pieces every branch needs, so that a consumer can
read a penalty's size, its hyperparameter names and their bounds without
knowing which branch it holds.

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
  [`print()`](https://rdrr.io/r/base/print.html) and by consumers that
  report which penalty a block carries.

- map:

  The matrix \\D\\, of \\m\\ rows and \\q\\ columns, or `NULL` for the
  identity. A Matrix object is kept in its own storage; a diagonal map
  is what standardization comes to and is recognized by its class.

- n_coef:

  The number of coefficients \\q\\. A single whole number.

- params:

  The hyperparameter names, in the order every derivative list is keyed
  by. `character(0)` for a penalty with none.

- params_bounds:

  A named list, one entry per hyperparameter, each a numeric pair giving
  an **open** interval. A value at either endpoint is rejected, so
  `(0, Inf)` excludes zero.

- link_params:

  A named list, one linkfunctions7 link per hyperparameter, carrying
  that hyperparameter's own interval onto the whole real line.

- params_smooth:

  A logical vector, one entry per hyperparameter, `TRUE` where the value
  is differentiable in it.

## Value

An S7 object of class `penalty` carrying the seven properties above. The
class is abstract: an object of exactly this class answers
[`print()`](https://rdrr.io/r/base/print.html), and every generic that
computes something rejects it.

## What the properties mean

With \\q\\ coefficients and a map of \\m\\ rows, `beta` is a vector of
length \\q\\, \\D\beta\\ has length \\m\\, and every derivative in the
coefficients comes back at length \\q\\ or shape \\q \times q\\. A
`NULL` map is the identity, and then \\m = q\\ and no arithmetic is
done.

Each hyperparameter carries a linkfunctions7 link in `link_params`,
mapping its own open interval onto the whole real line. A penalty is
therefore optimizable on the unconstrained scale: a caller works in
\\\eta = g(\theta)\\, never has to police the bounds, and
[`penalty_grad_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
and its second-order siblings answer on either scale through their
`scale` argument.

`params_smooth` records which hyperparameters the value is
differentiable in. It is `TRUE` for every hyperparameter of every
shipped branch; the slot exists so that a branch with a
non-differentiable hyperparameter can say so, as distributions7 does for
a Laplace location.

## What a subclass owes

Construct this class directly only to write a branch of your own. The
generics registered on `penalty` itself are the refusals and the
defaults:
[`penalty_value()`](https://statmodels7.github.io/penalties7/reference/penalty_value.md)
has no method here at all and a bare `penalty` object rejects it, while
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md),
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md),
[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
and the three quantities beside it answer `FALSE` or reject. A subclass
supplies the value, the gradient, the Hessian, the three hyperparameter
blocks and
[`penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.md),
and
[`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
then says whether they agree with each other.

The four branches that ship are
[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md),
[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md),
[`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md)
with
[`mcp_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md),
and
[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md),
with
[`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md)
a fifth built on a parameters7 matrix parameter.

## See also

[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md),
[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md),
[`scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.md),
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
# Every branch inherits from this class, so a consumer can test for it.
pen <- quadratic_penalty(crossprod(diff(diag(4))), map = NULL)
S7::S7_inherits(pen, penalty)
#> [1] TRUE

# The properties a consumer reads without knowing the branch.
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

# A second-difference penalty over four coefficients, restricted to the
# first three by a map: three coefficients in, two rows out.
D <- diff(diag(3))
mapped <- quadratic_penalty(diag(2), map = D)
mapped
#> quadratic penalty on 3 coefficient(s) through 2 row(s); theta: lambda
```
