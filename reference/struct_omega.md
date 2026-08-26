# Precision, and Its Derivatives, From the Structure

The four helpers that read the prior's precision and its derivatives off
the structure, transporting from the covariance where that is what the
structure describes. `struct_omega()` returns \\\Omega\\, `struct_d1()`
its first derivatives, `struct_d2()` its second, and `struct_logdet()`
the log-determinant with as many orders as asked for.

## Usage

``` r
struct_omega(pen, eta)

struct_d1(pen, eta, omega = NULL)

struct_d2(pen, eta, omega = NULL)

struct_logdet(pen, eta, order = 2L)
```

## Arguments

- pen:

  A
  [`StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/StructuredPenalty.md)
  object.

- eta:

  The structure's free vector, as
  [`struct_eta()`](https://statmodels7.github.io/penalties7/reference/struct_eta.md)
  returns it.

- omega:

  The precision, when the caller already has it, so that a covariance
  structure is not inverted twice. `NULL` computes it. `struct_d1()` and
  `struct_d2()` only.

- order:

  The highest derivative wanted, 0, 1 or 2.

## Value

`struct_omega()` a symmetric base matrix of side
`pen@structure@dimension`. `struct_d1()` a list of such matrices, one
per free value, in `pen@params` order. `struct_d2()` a list of such
matrices, one per unordered pair, keyed by the colon-joined sorted
names. `struct_logdet()` a list with `value` (a single number) and,
according to `order`, `d1` (a named numeric vector, one per free value)
and `d2` (a named list, one number per pair, keyed as `struct_d2()` is).

## Details

The branch is written in the precision throughout, so the transport
happens here and every method is the same arithmetic in both roles.
Writing it twice would be two implementations of one prior.

For a covariance structure, with \\A_k\\ and \\A\_{kl}\\ the structure's
own derivative arrays,

\$\$\Omega = \Sigma^{-}, \qquad \partial_k\Omega = -\Omega A_k \Omega,
\qquad \partial\_{kl}\Omega = \Omega\left(A_k\Omega A_l + A_l\Omega
A_k\right)\Omega - \Omega A\_{kl}\Omega,\$\$

and the log-determinant is negated at every order. For a precision
structure each helper unwraps the structure's answer and returns it.

Second-order components are keyed as parameters7 keys them: the two free
names joined by a colon, the pair sorted by position in `free_names`.
The key is constructed from the pair, never parsed out of a name, since
a free name may itself contain a colon.

## See also

[`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md)
for the prior these serve,
[`parameters7::param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.html)
and
[`parameters7::param_dlogdet()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.html)
for the contract they read.
