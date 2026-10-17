# Contracting a Multivariate Parent's Response Tensor Per Block

`dp_contract()` contracts the trailing response indices of a per-block
tensor against one direction each: given an array of dimension \\p
\times \dots \times p \times n\\ and \\m\\ matrices of \\n\\ rows and
\\p\\ columns, it returns the \\p \times p \times n\\ array whose slice
\\i\\ is \\\sum\_{c_1 \dots c_m} A\_{ab c_1 \dots c_m
i}\\v\_{1,ic_1}\cdots v\_{m,ic_m}\\. That is the derivative of the
response Hessian of block \\i\\ along the coordinates of the block in
the directions, which
[`dp_blockdiag()`](https://statmodels7.github.io/penalties7/reference/dp_blockdiag.md)
then places on the diagonal. The function handles arrays of four and
five dimensions, that is one or two directions.

`dp_parent_or_reject()` calls the generic `fun` of the parent. If that
call signals an error of any kind, it signals in its place an error that
names the penalty and the derivative, so that a criterion receives
either the complete tensor or an error.

## Usage

``` r
dp_contract(A, dirs)

dp_parent_or_reject(pen, what, fun, a, theta)
```

## Arguments

- A:

  A numeric array with \\n\\ as its last dimension.

- dirs:

  A list of \\n \times p\\ matrices, one per contracted index; the last
  remaining index of the array is contracted first, against the last
  direction.

- pen:

  A
  [`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
  object.

- what:

  A word naming the derivative, for the message.

- fun:

  The distributions7 generic to call.

- a:

  The argument of the parent, from
  [`dp_arg()`](https://statmodels7.github.io/penalties7/reference/dp_arg.md).

- theta:

  A named list of the free parameters of the parent.

## Value

`dp_contract()` returns a \\p \times p \times n\\ numeric array.
`dp_parent_or_reject()` returns the value that `fun` returns, or signals
an error.

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
