# Precision, and Its Derivatives, From the Structure

This page documents the four helpers that read the prior's precision and
its derivatives off the structure. `struct_omega()` returns \\\Omega\\,
`struct_d1()` its first derivatives, `struct_d2()` its second, and
`struct_logdet()` the log-determinant with its derivatives up to the
requested order.

## Usage

``` r
struct_omega(penalty, eta)

struct_d1(penalty, eta, omega = NULL)

struct_d2(penalty, eta, omega = NULL)

struct_logdet(penalty, eta, order = 2L)
```

## Arguments

- penalty:

  A
  [`StructuredPenalty()`](https://statmodels7.github.io/penalties7/reference/StructuredPenalty.md)
  object.

- eta:

  The structure's free vector, as
  [`struct_eta()`](https://statmodels7.github.io/penalties7/reference/struct_eta.md)
  returns it.

- omega:

  Accepted and unused; the calls in this package never pass it. Only
  `struct_d1()` and `struct_d2()` take it.

- order:

  The highest log-determinant derivative wanted: `0`, `1` or `2`. Only
  `struct_logdet()` takes it, and its default is `2L`.

## Value

`struct_omega()` returns a symmetric base matrix of side
`penalty@structure@dimension`. `struct_d1()` returns a list of such
matrices, one per free value, in `penalty@params` order. `struct_d2()`
returns a list of such matrices, one per unordered pair, keyed by the
colon-joined sorted names. `struct_logdet()` returns a list with `value`
(a single number) and, according to `order`, `d1` (a named numeric
vector, one per free value) and `d2` (a named list, one number per pair,
keyed as `struct_d2()` is).

## Details

The structure is the precision, so each helper removes the class from
the structure's own result, returns that result and does not convert
between covariance and precision. To read a structure `s` as a
covariance, pass `parameters7::inverse_of(s)`, which carries
\\\partial_k\Omega = -\Omega B_k \Omega\\, with \\B_k\\ the derivative
of the covariance, and its higher orders inside the structure, where the
log-determinant stays exact.

Second-order components are keyed as parameters7 keys them: a key joins
the two free names with a colon, and the pair is sorted by position in
`free_names`. The key is constructed from the pair, never parsed out of
a name, since a free name may itself contain a colon.

## See also

[`structured_penalty()`](https://statmodels7.github.io/penalties7/reference/structured_penalty.md)
for the prior that these serve,
[`parameters7::param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.html)
and
[`parameters7::param_dlogdet()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.html)
for the generics that they call.
