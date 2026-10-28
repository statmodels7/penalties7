# The Value of a Parent's Parameter, Free or Held

Returns the parent distribution's parameter `name`: from the aligned
hyperparameter list when the penalty carries it as a free
hyperparameter, and from the parent's `fixed_params` when
[`distributions7::fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.html)
holds it at a value.

## Usage

``` r
.prox_param(penalty, theta, name)
```

## Arguments

- penalty:

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

The closed forms need the value of the parameter, wherever it is stored.
A penalty whose scale is held by
[`distributions7::fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.html)
is the same penalty as one whose scale is free, evaluated at that scale,
and
[`penalty_value()`](https://statmodels7.github.io/penalties7/reference/penalty_value.md),
[`penalty_gradient()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
and
[`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
return the same results for the two.
[`penalty_prox()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.md)
therefore reads the value from `theta` when it is a free hyperparameter
and from the parent's `fixed_params` otherwise.

## Errors

The message names the penalty and the parameter when the parameter is
neither a free hyperparameter nor held by
[`distributions7::fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.html).

## See also

[`penalty_prox.DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_prox.DistribPenalty.md),
[`penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.md)
