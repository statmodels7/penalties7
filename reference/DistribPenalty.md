# S7 Class for the Separable Penalty

The class
[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md)
builds: a penalty formed by applying a distributions7 log-density to the
successive blocks of \\D\beta\\. Beyond the seven properties every
penalty carries it stores the parent distribution, the points where the
parent's log-density has a kink, and the block width, which is one for a
univariate parent and the dimension for a multivariate one.

## Usage

``` r
DistribPenalty(
  penalty_name = character(0),
  map = NULL,
  n_coef = integer(0),
  params = character(0),
  params_bounds = list(),
  link_params = list(),
  params_smooth = logical(0),
  parent = NULL,
  kinks = integer(0),
  block = integer(0)
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

- parent:

  The distributions7 object the log-density comes from, usually a
  [`distributions7::fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.html)
  wrapper holding the location at zero.

- kinks:

  A numeric vector of the points where the parent's log-density is not
  differentiable in its argument, possibly empty. Derived by
  [`distrib_kinks()`](https://statmodels7.github.io/penalties7/reference/distrib_kinks.md)
  unless the constructor was told otherwise.

- block:

  The block width: `1L` for a univariate parent and the parent's `n_dim`
  for a multivariate one.

## Value

An S7 object of class `DistribPenalty`, inheriting from
[`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md),
with the seven inherited properties and the three above.

## Details

The hyperparameters are the parent's own free parameters, and `params`,
`params_bounds` and `link_params` are read straight off the distribution
object. A caller who knows the parent therefore knows the penalty's
hyperparameters, their bounds and their links without being told again.

## See also

[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md)
for the constructor to use,
[`ridge_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.md)
and its siblings for the named instances,
[`penalty_value.DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.DistribPenalty.md)
for what the branch computes.

## Examples

``` r
pen <- lasso_penalty(n_coef = 2)
S7::S7_inherits(pen, DistribPenalty)
#> [1] TRUE

# The parent is a fixed Laplace, the kink is at zero, blocks are single
# coordinates, and the hyperparameter is the parent's own rate.
pen@parent
#> Distribution: Fixed Laplace2 [mu=0]
#> Type:         Continuous
#> Dimensions:   univariate
#> 
#> Parameters:
#>   lambda (rate)               | Link: log        | Domain: (0, Inf)
#> 
#> Fixed:
#>   mu = 0
pen@kinks
#> [1] 0
pen@block
#> [1] 1
pen@params
#> [1] "lambda"
```
