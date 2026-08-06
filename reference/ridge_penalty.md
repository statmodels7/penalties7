# Named Separable Penalties

The canonical instances of
[`distrib_penalty`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md),
shipped as constructors so the model layer can name what it means. Ridge
is the Gaussian at zero with the scale free; the lasso is the Laplace at
zero, with its kink declared; the heavy-tailed prior is the Student t at
zero, whose \\\nu\\ is estimable exactly because the normalizing
constant is kept.

## Usage

``` r
ridge_penalty(map = NULL, n_coef = 1L)

lasso_penalty(map = NULL, n_coef = 1L)

heavy_penalty(map = NULL, n_coef = 1L)
```

## Arguments

- map:

  The matrix \\D\\, or `NULL` (default) for the identity.

- n_coef:

  The number of coefficients when `map` is `NULL`.

## Value

An object of class `DistribPenalty`.

## Examples

``` r
pen <- ridge_penalty(n_coef = 2)
penalty_gradient(pen, c(1, -1), list(sigma = 1))
#> [1]  1 -1
```
