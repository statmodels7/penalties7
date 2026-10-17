# Carry a Parent's Paired Components Into Coefficient Space

Re-keys a distributions7 component keyed by parameter pair into this
package's own keys, and places each through the penalty's map.

## Usage

``` r
carry_pairs(pen, comp, carry)
```

## Arguments

- pen:

  A
  [`penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- comp:

  The components of the parent, keyed by parameter pair.

- carry:

  A function that places one component into coefficient space.

## Value

A named list keyed by hyperparameter pair, in the order of
[`ptheta_pairs()`](https://statmodels7.github.io/penalties7/reference/ptheta_pairs.md).

## Details

The two enumerations of pairs are built the same way from the same
names, so a key of one is a key of the other. Each pair is looked up by
name, first in `pen@params` order and then in the reversed order, and
not by position or by splitting a key, since a hyperparameter whose own
name contains the separator could not be split back into its two names.
If the pair is absent in both orders, the function signals an error that
names the pair.
