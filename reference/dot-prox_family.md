# The Family Name Under a Separable Penalty's Wrappers

Returns the bare name of the distribution on which a separable penalty
is built, with a `fixed()` wrapper and any bracketed parameter list
removed, so that the closed-form branches can be recognized by the
family underneath.

## Usage

``` r
.prox_family(pen)
```

## Arguments

- pen:

  A
  [`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
  object.

## Value

A single string: `"gaussian1"`, `"laplace2"`, `"laplace"` or `"enet"`
for the families with a closed form (the Laplace in two
parametrizations), and the name of the parent's family otherwise.

## Details

The parent of a separable penalty is usually wrapped: the lasso's parent
is `fixed laplace2 [mu=0]`, and a separable Gaussian prior built as
`distrib_penalty(fixed(gaussian1_distrib(), mu = 0), ...)` has the
parent `fixed gaussian1 [mu=0]`. A match on the raw `distrib_name` would
miss both. Removing the leading `fixed ` and everything from the first
bracket leaves `laplace2` and `gaussian1`.

The name is used only to reach a closed form. Whether the parent is
centered at zero is tested separately, by evaluating the penalty's
gradient at the origin, because the family name does not carry the
location.

## See also

[`penalty_prox.DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.DistribPenalty.md),
[`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
