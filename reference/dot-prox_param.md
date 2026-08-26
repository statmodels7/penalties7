# The Value of a Parent's Parameter, Free or Held

Returns the parent distribution's parameter `name`: from the aligned
hyperparameter list when the penalty carries it as a free
hyperparameter, and from the parent's `fixed_params` when
[`distributions7::fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.html)
holds it at a value.

## Usage

``` r
.prox_param(pen, theta, name)
```

## Arguments

- pen:

  A
  [`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
  object.

- theta:

  The aligned hyperparameter list, as
  [`align_ptheta()`](https://statmodels7.github.io/penalties7/reference/align_ptheta.md)
  returns it.

- name:

  The parameter's name, a single string.

## Value

A single number.

## Details

The closed forms below need the number, not the place it is kept. A
penalty whose scale is held is the same penalty as one whose scale is
free, read at that scale:
[`penalty_value()`](https://statmodels7.github.io/penalties7/reference/penalty_value.md),
[`penalty_gradient()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
and
[`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
all already answer identically for the two. Reading only `theta` made
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
the one generic that did not, and it stopped with
`attempt to select less than one element in get1index` rather than
saying so, because [`which()`](https://rdrr.io/r/base/which.html) of an
empty comparison is `integer(0)`.

## See also

[`penalty_prox.DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.DistribPenalty.md),
[`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
