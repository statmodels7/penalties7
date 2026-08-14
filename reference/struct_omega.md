# Precision, and Its Derivatives, From the Structure

The prior's precision and its first and second derivatives in the
structure's free values, transported from the covariance where that is
what the structure describes.

## Usage

``` r
struct_omega(pen, eta)

struct_d1(pen, eta, omega = NULL)

struct_d2(pen, eta, omega = NULL)

struct_logdet(pen, eta, order = 2L)
```

## Arguments

- pen:

  A `StructuredPenalty` object.

- eta:

  The structure's free vector.

- omega:

  The precision, when the caller already has it.

- order:

  The highest derivative wanted, 0, 1 or 2.

## Value

A matrix (`struct_omega`) or a list of matrices.
