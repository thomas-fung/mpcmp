# Functions to Compute Various Expected Values for the COM-Poisson Distribution

Functions to approximate the various expected values for the COM-Poisson
distribution via truncation. The standard COM-Poisson parametrization is
being used here. The lambda and nu values are recycled to match the
length of the longer one and that would determine the length of the
results.

## Usage

``` r
comp_mean_logfactorialy(lambda, nu, log.Z, summax = 100)

comp_mean_ylogfactorialy(lambda, nu, log.Z, summax = 100)

comp_means(lambda, nu, log.Z, summax = 100)

comp_variances(lambda, nu, log.Z, summax = 100)

comp_variances_logfactorialy(lambda, nu, log.Z, summax = 100)
```

## Arguments

- lambda, nu, :

  rate and dispersion parameters. Must be positives.

- log.Z, :

  an optional vector specifying normalizing constant Z in log scale.

- summax:

  maximum number of terms to be considered in the truncated sum. The
  default is to sum to 100.

## Value

`comp_mean_logfactorialy` gives the mean of *log(Y!)*.

`comp_mean_ylogfactorialy` gives the mean of *ylog(Y!)*.

`comp_means` gives the mean of *Y*.

`comp_variances` gives the variance of *Y*.

`comp_variances_logfactorialy` gives the variance of *log(Y!)*.
