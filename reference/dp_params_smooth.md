# Helpers of the Separable Constructor

`dp_params_smooth()` returns the parent's `params_smooth` aligned to its
`params`, with `TRUE` for a parameter that `params_smooth` does not
list. `dp_is_multinomial()` returns `TRUE` when the object, or any
object that it wraps through `parent_distrib`, is a distributions7
multinomial.

## Usage

``` r
dp_params_smooth(d)

dp_is_multinomial(d)
```

## Arguments

- d:

  A distributions7 object.

## Value

`dp_params_smooth()` a named logical vector, one entry per element of
`d@params`. `dp_is_multinomial()` a single logical.

## See also

[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md)
