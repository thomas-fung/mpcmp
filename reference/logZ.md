# Calculate the Normalizing Constant in log scale for COM-Poisson distribution

A function to approximate the normalizing constant for COM-Poisson
distributions via truncation. The standard COM-Poisson parametrization
is being used here.

## Usage

``` r
logZ(log_lambda, nu, summax = 100)
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
normalizing constant \\\log Z(\lambda, \nu)\\, recycled to the length of
the longer of `log_lambda` and `nu`.

## Details

As of version 0.2.0 of this package, `logZ` will supersede `Z` for
calculating the normalizing constant. `logZ` utilised a method that can
calculate `log(exp(logx) + exp(logy))` in a somewhat numerically stable
way.

This function was originally purposed in the `cmpreg` package of Ribeiro
Jr, Zeviani & Demétrio (2019).

## References

Ribeiro Jr, E. E., Zeviani, W. M., Demétrio, C. G. B. (2019) `cmpreg`:
Reparametrized COM-Poisson Regression Models. R package version 0.0.1.
