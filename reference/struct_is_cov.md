# Whether the Structure Describes the Covariance

Reads the structure's declared role and answers `TRUE` for
`"covariance"`. Every method of the branch is written in the precision,
so this is the one place that decides whether a transport is needed.

## Usage

``` r
struct_is_cov(pen)
```

## Arguments

- pen:

  A
  [`StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/StructuredPenalty.md)
  object.

## Value

A single logical. `FALSE` for a structure of role `"precision"`, which
is the only other value
[`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md)
admits.

## See also

[`struct_omega()`](https://statmodels7.github.io/penalties7/reference/struct_omega.md)
for the transport this gates
