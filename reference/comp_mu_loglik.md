# Calculate the Log-Likelihood of the COM-Poisson model

A function to compute the log-likelihood of the COM-Poisson model.

## Usage

``` r
comp_mu_loglik(param, y, xx, offset, summax)

comp_mu_neg_loglik_log_nu_only(log_nu, mu, y, summax)
```

## Arguments

- param:

  numeric vector: the model coefficients & the current value of `nu`. It
  is assumed that `nu` is in the last position of `param`.

- y:

  numeric vector: response variable

- xx:

  numeric matrix: the explanatory variables

- offset:

  numeric vector: a vector of length equal to the number of cases

- summax:

  maximum number of terms to be considered in the truncated sum

- log_nu:

  numeric: nu in log-scale

- mu:

  numeric vector: fitted mean parameters

## Value

`comp_mu_loglik` returns the log-likelihood value of the COM-Poisson
model based on Huang (2017). `comp_mu_neg_loglik_log_nu_only` returns
the negative log-likelihood value of the COM-Poisson model based on
Ribeiro Jr et al. (2018)'s specification to use in conjunction with
`optim`.
