# Does a Penalty Supply a Proximal Operator?

`TRUE` when
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
can be evaluated for this penalty, so that a caller may route a block to
a proximal method without provoking an error. `FALSE` when the penalty
has none, and then the block belongs to a smooth method or to a scheme
of its own.

## Usage

``` r
has_prox(pen)
```

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object of any branch. Anything else is rejected.

## Value

A single logical.

## What it tests

The predicate asks four things in turn: whether the class registers a
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
method of its own rather than inheriting the base class's refusal;
whether the penalty is quadratic or structured, in which case the
operator is a linear solve at any map; whether a separable parent is
read in blocks, whose coordinates do not separate; and otherwise whether
the map is absent or diagonal.

Measured over the shipped branches:

|                                         |              |
|-----------------------------------------|--------------|
| penalty                                 | `has_prox()` |
| quadratic, at any map                   | `TRUE`       |
| structured                              | `TRUE`       |
| ridge, lasso, elastic net, heavy-tailed | `TRUE`       |
| lasso under a **diagonal** map          | `TRUE`       |
| SCAD, MCP                               | `TRUE`       |
| **additive**                            | `FALSE`      |
| lasso under a **general** map           | `FALSE`      |
| a separable penalty read in blocks      | `FALSE`      |

## What `TRUE` does not promise

It says the penalty has an operator, not that every step reaches it.
SCAD and MCP answer `TRUE` and still reject a step at or beyond \\a -
1\\ and \\\gamma\\, where the subproblem is not convex and the operator
is set-valued. A Laplace or elastic-net parent that is not centered at
zero answers `TRUE` and rejects at the call. The predicate is about the
penalty; those two conditions are about the step and the parent.

## See also

[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
for the operator,
[`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
for it as a table,
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
and
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
for the other two predicates a consumer routes on.

## Examples

``` r
c(lasso = has_prox(lasso_penalty()),
  scad  = has_prox(scad_penalty()),
  ridge = has_prox(ridge_penalty()))
#> lasso  scad ridge 
#>  TRUE  TRUE  TRUE 

# The additive branch is the one that ships without an operator.
has_prox(additive_penalty(list(diag(3), diag(c(1, 1, 0)))))
#> [1] FALSE

# A map decides it for a separable penalty: diagonal yes, general no.
c(diagonal = has_prox(lasso_penalty(map = Matrix::Diagonal(x = c(1, 2, 3)))),
  general  = has_prox(lasso_penalty(map = rbind(c(1, -1, 0), c(0, 1, -1)))))
#> diagonal  general 
#>     TRUE    FALSE 

# TRUE does not mean every step works: SCAD needs step < a - 1.
has_prox(scad_penalty(n_coef = 1))
#> [1] TRUE
try(penalty_prox(scad_penalty(n_coef = 1), 2, 2.7, list(lambda = 1, a = 3.7)))
#> Error : the SCAD proximal operator needs step < a - 1 (2.7 here): beyond that the
#>   subproblem is not convex and the operator is set-valued.
```
