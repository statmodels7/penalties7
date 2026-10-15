# Drawing From an Additive Penalty

One draw from the Gaussian prior with precision \\S(\lambda)\\, block by
block where the penalty is repeated blockwise.

## Arguments

- pen:

  An
  [`AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/AdditivePenalty.md)
  object.

- theta:

  A named list of `lambda1`, `lambda2`, ...

- ...:

  Unused.

## Value

A numeric vector of length `pen@n_coef`.

## Details

With one block the draw is the base method's. Repeated blockwise, the
precision is \\I_m \otimes S_1(\lambda)\\, so one triangular factor of
\\S_1\\ gives every block: with \\S_1 = R^\top R\\ the vector
\\R^{-1}z\\ has covariance \\S_1^{-1}\\. A block whose precision is
singular is not a proper prior and comes back `NA`, as the base method
returns it.

## See also

[`penalty_draw()`](https://statmodels7.github.io/penalties7/reference/penalty_draw.md)

## Examples

``` r
pen <- additive_penalty(list(diag(2), diag(c(1, 2))), blocks = 3)
set.seed(1)
penalty_draw(pen, list(lambda1 = 1, lambda2 = 4))
#> [1] -0.28015866  0.06121444 -0.37370448  0.53176027  0.14736036 -0.27348946
```
