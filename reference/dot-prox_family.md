# The Family Name Under a Separable Penalty's Wrappers

Returns the bare name of the distribution a separable penalty is built
on, with a `fixed()` wrapper and any bracketed parameter list stripped
off, so that the closed-form branches can be recognized by the family
underneath.

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

A single string: `"gaussian1"`, `"laplace2"`, `"laplace"`, `"enet"` for
the four closed-form families, and whatever the parent is called
otherwise.

## Details

A separable penalty's parent is typically wrapped: the lasso is
`fixed laplace2 [mu=0]` and the ridge's separable twin is
`fixed gaussian1 [mu=0]`. Matching on the raw `distrib_name` would miss
both. Stripping the leading `fixed ` and everything from the first
bracket leaves `laplace2` and `gaussian1`.

The name is used only to reach a closed form. Whether the parent is
centered where the quadratic pull is is asked separately, by evaluating
the gradient at the origin, since a family name does not say where a
location sits.

## See also

[`penalty_prox.DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.DistribPenalty.md),
[`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
