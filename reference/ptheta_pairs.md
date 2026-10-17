# The Hyperparameter Pair Names

Returns the keys of a penalty's second hyperparameter derivatives and
the pairs that the keys represent. The \\p\\ diagonal keys come first,
in `params` order, the \\p(p-1)/2\\ upper off-diagonal pairs follow, and
each key joins the two names with an underscore. For
`c("lambda", "alpha")` the keys are `lambda_lambda`, `alpha_alpha`,
`lambda_alpha`. This is the order that
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
returns and that
[`ptheta_to_link()`](https://statmodels7.github.io/penalties7/reference/ptheta_to_link.md)
reads.

## Usage

``` r
ptheta_pairs(params)
```

## Arguments

- params:

  A character vector of hyperparameter names, in the order in which the
  penalty holds them. The vector may be empty.

## Value

A named list of length \\p(p+1)/2\\, empty when `params` is. Each
element is a character pair naming the two hyperparameters of the
derivative, and each name joins those two with an underscore.

## Details

The diagonals come first, rather than in lexicographic order, so that
the first \\p\\ entries are the diagonal of the Hessian. The same
convention names the Hessian components of distributions7.

`character(0)` gives an empty named list. A penalty with no free
hyperparameters has no pair in which to differentiate, and the named
empty list is the shape that
[`penalty_grad_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
and
[`penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
return in that case. The function treats this case separately, because
`paste0(character(0), "_", character(0))` returns the single string
`"_"` and [`stats::setNames()`](https://rdrr.io/r/stats/setNames.html)
would then signal an error on the empty list of pairs.
