# A Point to Read the Parent At

Returns one value per free parameter, at which the parent can be
evaluated when no hyperparameters are in force. Used by
[`distrib_kinks()`](https://statmodels7.github.io/penalties7/reference/distrib_kinks.md),
which has to differentiate the log-density before any penalty exists.

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

The rule is the toolkit's probe rule: the midpoint of a parameter's
bounds where both are finite, one above the lower bound where only that
is finite, one below the upper where only that is, and `1` where neither
is. So a scale on \\(0, \infty)\\ is probed at 1 and a mixing weight on
\\(0, 1)\\ at 0.5.

Whether a point is a kink does not depend on where the other parameters
sit, so any admissible point serves and the rule only has to give one.

## See also

[`distrib_kinks()`](https://statmodels7.github.io/penalties7/reference/distrib_kinks.md),
[`has_jump()`](https://statmodels7.github.io/penalties7/reference/has_jump.md)
