# Extract COM-Poisson Model Residuals

`residuals` is a generic function which extracts model residuals from
objects returned by the modelling function `glm.comp`. `resid` is an
alias for `residuals` .

## Usage

``` r
# S3 method for class 'cmp'
residuals(object, type = c("deviance", "pearson", "response"), ...)
```

## Arguments

- object:

  an object class 'cmp', obtained from a call to `glm.cmp`.

- type:

  the `type` of residuals which should be returned. The alternatives
  are: 'deviance' (default), 'pearson' and 'response'. Can be
  abbreviated.

- ...:

  other arguments passed to or from other methods (currently unused).

## Value

Residuals extracted from the object `object`.

## See also

[`coef.cmp`](https://thomas-fung.github.io/mpcmp/reference/coef.cmp.md),
[`fitted.cmp`](https://thomas-fung.github.io/mpcmp/reference/fitted.cmp.md),
[`glm.cmp`](https://thomas-fung.github.io/mpcmp/reference/glm.cmp.md)
