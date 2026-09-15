# Test for a whole number

Test for integer/whole number vector

## Usage

``` r
is.wholenumber(x, tol = .Machine$double.eps^0.5)
```

## Arguments

- x:

  numeric vector to be tested

- tol:

  numeric; precision level

## Value

A logical vector, the same length as `x`, indicating whether each
element is (within tolerance `tol` of) a whole number.
