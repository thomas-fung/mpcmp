# Calculate the Normalizing Constant in log scale for COM-Poisson distribution The calculation of the function `logZ` will be performed here. This function is used to approximate the normalizing constant for COM-Poisson distributions via truncation. The standard COM-Poisson parametrization is being used here.

It is assumed that vectors log_lambda & nu are of equal length.

## Usage

``` r
logZ_c(log_lambda, nu, summax)
```

## Arguments

- log_lambda:

  rate parameter in log scale.

- nu:

  dispersion parameter, straightly positive.

- summax:

  maximum number of terms to be considered in the truncated sum.

## Value

A numeric vector giving the (approximate) log of the COM-Poisson
normalizing constant \\\log Z(\lambda, \nu)\\ for each element of
`log_lambda` and `nu`.

## Details

This function was originally purposed in the `cmpreg` package of Ribeiro
Jr, Zeviani & Demétrio (2019).

## References

Ribeiro Jr, E. E., Zeviani, W. M., Demétrio, C. G. B. (2019) `cmpreg`:
Reparametrized COM-Poisson Regression Models. R package version 0.0.1.
