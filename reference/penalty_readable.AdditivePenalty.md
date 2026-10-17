# Readable Quantities of an Additive Penalty

Returns `NULL` for a penalty built with `readable = "lambda"`, whose
smoothing parameters are the quantities reported, and the standard
deviations \\\sigma_k = \lambda_k^{-1/2}\\ for one built with
`readable = "sd"`.

## Arguments

- pen:

  An
  [`AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/AdditivePenalty.md)
  object.

- theta:

  A named list of `lambda1`, `lambda2`, ...

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

`NULL`, or a list of `value` (named `sigma1`, `sigma2`, ...),
`jacobian`, `transform` and `block`, as
[`penalty_readable()`](https://statmodels7.github.io/penalties7/reference/penalty_readable.md)
documents.

## Details

Each standard deviation depends only on its own parameter, so the
Jacobian is diagonal, with entries \\-\tfrac{1}{2}\lambda_k^{-3/2}\\,
and each interval is built on the log scale.

## See also

[`penalty_readable()`](https://statmodels7.github.io/penalties7/reference/penalty_readable.md)
for the generic,
[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
for the argument.

## Examples

``` r
pen <- additive_penalty(list(diag(2), diag(c(1, 2))), readable = "sd")
penalty_readable(pen, list(lambda1 = 4, lambda2 = 0.25))
#> $value
#> sigma1 sigma2 
#>    0.5    2.0 
#> 
#> $jacobian
#>        lambda1 lambda2
#> sigma1 -0.0625       0
#> sigma2  0.0000      -4
#> 
#> $transform
#> [1] "log" "log"
#> 
#> $block
#> [1] "sd" "sd"
#> 
```
