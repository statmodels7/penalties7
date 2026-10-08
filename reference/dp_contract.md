# Contracting a Multivariate Parent's Response Tensor Per Block

`dp_contract()` contracts the trailing response indices of a per-block
tensor against one direction each: given an array of dimension \\p
\times \dots \times p \times n\\ and \\m\\ matrices of \\n\\ rows and
\\p\\ columns, it returns the \\p \times p \times n\\ array whose slice
\\i\\ is \\\sum\_{c_1 \dots c_m} A\_{ab c_1 \dots c_m
i}\\v\_{1,ic_1}\cdots v\_{m,ic_m}\\. That is the movement of a block's
response Hessian along the block's coordinates of the directions, which
[`dp_blockdiag()`](https://statmodels7.github.io/penalties7/reference/dp_blockdiag.md)
then places on the diagonal.

`dp_parent_or_reject()` asks the parent for a response tensor and, where
the parent does not supply one, rejects with the penalty named, so that
a criterion is left without an exact gradient rather than given one
missing a piece.

## Usage

``` r
dp_contract(A, dirs)

dp_parent_or_reject(pen, what, fun, a, theta)
```

## Arguments

- A:

  A numeric array with \\n\\ as its last dimension.

- dirs:

  A list of \\n \times p\\ matrices, one per contracted index, the last
  index first contracted against the last direction.

- pen:

  A
  [`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
  object.

- what:

  A word naming the derivative, for the message.

- fun:

  The distributions7 generic to call.

- a:

  The parent's argument, from
  [`dp_arg()`](https://statmodels7.github.io/penalties7/reference/dp_arg.md).

- theta:

  A named list of the parent's free parameters.

## Value

`dp_contract()` a \\p \times p \times n\\ numeric array.
`dp_parent_or_reject()` the parent's answer, or an error.

## See also

[`penalty_dhessian_beta()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian_beta.md),
[`penalty_d2hessian_beta()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.md).

## Examples

``` r
A <- array(1, c(2, 2, 2, 3))
v <- matrix(c(1, 0, 2, 1, 1, 0), 3, 2)
dim(penalties7:::dp_contract(A, list(v)))
#> [1] 2 2 3
```
