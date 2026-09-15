# Extract the (Maximized) Log-Likelihood from a COM-Poisson Model Fit

An accessor function used to extract the (maximized) log-likelihood from
a 'cmp' object.

## Usage

``` r
# S3 method for class 'cmp'
logLik(object, ...)

# S3 method for class 'logLik.cmp'
print(x, ...)
```

## Arguments

- object:

  an object of class 'cmp' object, obtained from a call to `glm.cmp`

- ...:

  other arguments passed to or from other methods (currently unused).

- x:

  an object of class 'logLik.cmp', obtained from a call to `logLik.cmp`.

## Value

`logLik.cmp` returns an object of class `"logLik.cmp"`: the maximized
log-likelihood value with a `"df"` attribute giving the number of
estimated parameters, analogous to
[`logLik`](https://rdrr.io/r/stats/logLik.html).

`print.logLik.cmp` is called for its side effect of printing and returns
`x` invisibly.

## See also

[`coef.cmp`](https://thomas-fung.github.io/mpcmp/reference/coef.cmp.md),
[`fitted.cmp`](https://thomas-fung.github.io/mpcmp/reference/fitted.cmp.md),
[`glm.cmp`](https://thomas-fung.github.io/mpcmp/reference/glm.cmp.md)
