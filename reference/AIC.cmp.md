# Akaike's Information Criterion

A function calculating Akaike's Information Criterion (AIC) based on the
log-likelihood value extracted from
[`logLik.cmp`](https://thomas-fung.github.io/mpcmp/reference/logLik.cmp.md),
according to the formula *-2\\log-likelihood + k\\npar*, where *npar*
represents the number of parameters in the fitted model, and *k=2* for
the usual AIC or *k=log(n)* (*n* being the number of observations) for
the so-called BIC (Bayesian Information Criterion).

## Usage

``` r
# S3 method for class 'cmp'
AIC(object, ..., k = 2)
```

## Arguments

- object:

  an object class 'cmp' object, obtained from a call to `glm.cmp`

- ...:

  other arguments passed to or from other methods (currently unused).

- k:

  numeric: the *penalty* per parameter to be used; the default k = 2 is
  the classical AIC.

## Value

A numeric value with the corresponding AIC (or BIC, or ..., depends on
k).

## Details

When comparing models fitted by maximum likelihood to the same data, the
smaller the AIC or BIC, the better the fit.

## See also

[`logLik.cmp`](https://thomas-fung.github.io/mpcmp/reference/logLik.cmp.md),
[`nobs.cmp`](https://thomas-fung.github.io/mpcmp/reference/nobs.cmp.md),
[`glm.cmp`](https://thomas-fung.github.io/mpcmp/reference/glm.cmp.md)
