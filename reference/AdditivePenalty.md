# S7 Class for a Sum of Quadratic Penalties

The class
[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
builds: several quadratic penalties added together, each with a
smoothing parameter of its own. Beyond the seven properties every
penalty carries it stores the list of component matrices, already
carried through the map, and the rank and null basis of the sum.

## Usage

``` r
AdditivePenalty(
  penalty_name = character(0),
  map = NULL,
  n_coef = integer(0),
  params = character(0),
  params_bounds = list(),
  link_params = list(),
  params_smooth = logical(0),
  mats = list(),
  p_rank = integer(0),
  null_basis = NULL
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

- mats:

  A list of symmetric positive semidefinite matrices of the same side,
  already carried through the map, so each is `n_coef` by `n_coef`. One
  smoothing parameter goes with each, named `lambda1`, `lambda2` and so
  on in this order.

- p_rank:

  The rank of the sum, a single whole number, the same at every positive
  setting of the parameters.

- null_basis:

  An orthonormal basis of the components' shared null space, `n_coef` by
  `n_coef - p_rank`, read from the same eigendecomposition as `p_rank`
  and equally fixed.

## Value

An S7 object of class `AdditivePenalty`, inheriting from
[`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md),
with the seven inherited properties and the three above. Its `map` is
always `NULL`, a map given to the constructor having been absorbed into
`mats`.

## Details

The rank and the null basis are stored because both are properties of
the components alone: the null space of a sum of positive semidefinite
matrices is the intersection of the components' null spaces, so neither
moves with the parameters. Reading the rank off the assembled
\\S(\lambda)\\ instead would make it fall as the parameters spread
apart, which
[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
measures.

## See also

[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
for the constructor to use,
[`penalty_value.AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_value.AdditivePenalty.md)
for what the branch computes,
[`quadratic_penalty()`](https://statmodels7.github.io/penalties7/reference/quadratic_penalty.md)
for one matrix with one scale.

## Examples

``` r
pen <- additive_penalty(list(diag(3), diag(c(1, 0, 0))))
S7::S7_inherits(pen, AdditivePenalty)
#> [1] TRUE

# One hyperparameter per component, and a rank fixed at construction.
pen@params
#> [1] "lambda1" "lambda2"
c(rank = pen@p_rank, null_width = ncol(pen@null_basis))
#>       rank null_width 
#>          3          0 
length(pen@mats)
#> [1] 2
```
