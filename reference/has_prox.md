# Does a Penalty Supply a Proximal Operator?

Returns `TRUE` when
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
has an operator for this penalty, so that a caller can route a block to
a proximal method. Returns `FALSE` when the penalty has none, and the
block then goes to a smooth method or to a scheme of its own.

## Usage

``` r
has_prox(pen)
```

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object of any branch. Any other object signals an error.

## Value

A single logical.

## What it tests

The predicate tests four conditions in turn: whether the class registers
a
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
method of its own instead of inheriting the error method of the base
class; whether the penalty is quadratic or structured, in which case the
operator is a linear solve at any map; whether a separable parent is
read in blocks, in which case the coordinates do not separate and the
result is `FALSE`; and otherwise whether the map is absent or a Matrix
diagonal object with nonzero entries, as
[`map_diagonal()`](https://statmodels7.github.io/penalties7/reference/map_diagonal.md)
reads it.

The value for each shipped branch is:

|  |  |
|----|----|
| penalty | `has_prox()` |
| quadratic, at any map | `TRUE` |
| structured | `TRUE` |
| ridge, lasso, elastic net, heavy-tailed, with no map | `TRUE` |
| a separable penalty under a diagonal map with nonzero entries | `TRUE` |
| SCAD, MCP, with no map or a diagonal map with nonzero entries | `TRUE` |
| additive | `FALSE` |
| a separable, SCAD or MCP penalty under a general map | `FALSE` |
| a separable penalty read in blocks | `FALSE` |

## Conditions on the step and on the parent

A `TRUE` value means that the penalty has an operator; a call can still
signal an error for some steps and some parents. SCAD and MCP return
`TRUE` and reject a step for which \\t c_j\\ reaches \\a - 1\\ or
\\\gamma\\ in some coordinate (where \\c_j\\ is the curvature and \\t\\
is multiplied by \\d_j^2\\ under a diagonal map), because the subproblem
is then not convex and the operator is set-valued. A Laplace or
elastic-net parent that is not centered at zero, and a separable parent
that declares a kink and has no closed form, return `TRUE`, and the call
signals an error.

## See also

[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
for the operator,
[`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
for it as a table,
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
and
[`is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.md)
for the other two predicates that a fitting layer uses to route a block.

## Examples

``` r
c(lasso = has_prox(lasso_penalty()),
  scad  = has_prox(scad_penalty()),
  ridge = has_prox(ridge_penalty()))
#> lasso  scad ridge 
#>  TRUE  TRUE  TRUE 

# The additive branch has no operator.
has_prox(additive_penalty(list(diag(3), diag(c(1, 1, 0)))))
#> [1] FALSE

# For a separable penalty the value depends on the map: TRUE under a
# diagonal map, FALSE under a general one.
c(diagonal = has_prox(lasso_penalty(map = Matrix::Diagonal(x = c(1, 2, 3)))),
  general  = has_prox(lasso_penalty(map = rbind(c(1, -1, 0), c(0, 1, -1)))))
#> diagonal  general 
#>     TRUE    FALSE 

# SCAD has an operator, and a step at a - 1 signals an error.
has_prox(scad_penalty(n_coef = 1))
#> [1] TRUE
try(penalty_prox(scad_penalty(n_coef = 1), 2, 2.7, list(lambda = 1, a = 3.7)))
#> Error : the SCAD proximal operator needs step < a - 1 (2.7 here), where the
#>   step is multiplied by the curvature and, under a diagonal map, by
#>   the squared diagonal entry; the largest scaled step is 2.7. Beyond
#>   that bound the subproblem is not convex and the operator is
#>   set-valued.
```
