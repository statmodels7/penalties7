# The Hyperparameter Pair Names

The keys of a penalty's second hyperparameter derivatives, and the pairs
they stand for: the \\p\\ diagonals first, in `params` order, then the
\\p(p-1)/2\\ upper off-diagonal pairs, each joined by an underscore. For
`c("lambda", "alpha")` the keys are `lambda_lambda`, `alpha_alpha`,
`lambda_alpha`. This is the order
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
returns and the order
[`ptheta_to_link()`](https://statmodels7.github.io/penalties7/reference/ptheta_to_link.md)
reads.

## Usage

``` r
ptheta_pairs(params)
```

## Arguments

- params:

  A character vector of hyperparameter names, in the order the penalty
  holds them. May be empty.

## Value

A named list of length \\p(p+1)/2\\, empty when `params` is. Each
element is a character pair naming the two hyperparameters
differentiated in, and each name is those two joined by an underscore.

## Details

Diagonals first rather than lexicographically, because a consumer
reading only the variances can take the first \\p\\ entries. The same
convention names distributions7's Hessian components.

`character(0)` gives the empty **named** list. A penalty with no free
hyperparameters has no pair to differentiate in, and the named form
keeps
[`penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
the same shape as its two siblings
[`penalty_grad_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md)
and
[`penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.md),
both of which already answered for that case. The guard is needed rather
than incidental: `paste0(character(0), "_", character(0))` recycles the
zero-length argument against the length-one literal and gives the single
string `"_"`, so without it the names are one element long while the
list of pairs is empty and
[`stats::setNames()`](https://rdrr.io/r/stats/setNames.html) raises.
