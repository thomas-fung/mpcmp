# Calculate the Normalizing Constant for COM-Poisson distribution

A function to approximate the normalizing constant for COM-Poisson
distributions via truncation. The standard COM-Poisson parametrization
is being used here. Notice that the sum is hard coded to truncate at 100
so the approximation will be quite bad if the COM-Poisson has a large
rate or mean.

## Usage

``` r
Z(lambda, nu, log.z = FALSE, summax)
```

## Arguments

- lambda:

  rate parameter, straightly positive

- nu:

  dispersion parameter, straightly positive

- log.z:

  logical; if `TRUE`, normalising constant \\Z\\ are returned as
  \\log(Z)\\.

- summax:

  maximum number of terms to be considered in the truncated sum

## Value

A numeric vector giving the (approximate) COM-Poisson normalizing
constant \\Z(\lambda, \nu)\\ (or its natural log if `log.z = TRUE`),
recycled to the length of the longer of `lambda` and `nu`.
