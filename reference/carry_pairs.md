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
  [`penalty`](https://statmodels7.github.io/penalties7/reference/penalty.md)
  object.

- comp:

  The parent's components, keyed by parameter pair.

- carry:

  A function placing one component into coefficient space.

## Value

A named list keyed by hyperparameter pair.

## Details

The two enumerations of pairs are built the same way from the same
names, so a key from one is a key of the other; it is looked up by name
in both orders rather than by position, since a hyperparameter whose own
name contains the separator would not survive being taken apart.
