# Combine R Objects by Columns

Take a sequence of vector, matrix or data-frame arguments and combine
them by columns. `CBIND` is used within the package over `cbind` to
recycle the shorter arguments so that their number of rows would match.

## Usage

``` r
CBIND(..., deparse.level = 1)
```

## Arguments

- ...:

  (generalized) vectors or matrices. These can be given as named
  arguments

- deparse.level:

  integer; deparse.level = 0 constructs no labels, deparse.level = 1
  (the default) or \> 1 constructs labels from the arguments names.

## Value

A matrix formed by column-binding the (recycled) arguments; see
[`cbind`](https://rdrr.io/r/base/cbind.html).
