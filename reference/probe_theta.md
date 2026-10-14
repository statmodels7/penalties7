# A Point to Read the Parent At

Returns one value per free parameter, at which the parent can be
evaluated when no hyperparameters are given.
[`distrib_kinks()`](https://statmodels7.github.io/penalties7/reference/distrib_kinks.md)
uses it to differentiate the log-density before any penalty exists.

## Usage

``` r
probe_theta(d)
```

## Arguments

- d:

  A distributions7 object.

## Value

A named list of one number per element of `d@params`, in that order.
[`list()`](https://rdrr.io/r/base/list.html) for a distribution with no
free parameters.

## Details

The value of a parameter is the midpoint of its bounds where both are
finite, one above the lower bound where only that bound is finite, one
below the upper bound where only that bound is finite, and `1` where
neither is. A scale on \\(0, \infty)\\ is therefore probed at 1 and a
mixing weight on \\(0, 1)\\ at 0.5.

The rule assumes that whether a point is a kink does not depend on the
values of the other parameters, so that any admissible point serves.

## See also

[`distrib_kinks()`](https://statmodels7.github.io/penalties7/reference/distrib_kinks.md),
[`has_jump()`](https://statmodels7.github.io/penalties7/reference/has_jump.md)
