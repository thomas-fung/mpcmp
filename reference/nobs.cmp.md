# Extract the Number of Observation from a COM-Poisson Model Fit

An accessor function used to extract the number of observation from a
'cmp' object.

## Usage

``` r
# S3 method for class 'cmp'
nobs(object, ...)
```

## Arguments

- object:

  an object class 'cmp' object, obtained from a call to `glm.cmp`

- ...:

  other arguments passed to or from other methods (currently unused).

## Value

The number of observations extracted from the object `object`.

## See also

[`coef.cmp`](https://thomas-fung.github.io/mpcmp/reference/coef.cmp.md),
[`fitted.cmp`](https://thomas-fung.github.io/mpcmp/reference/fitted.cmp.md),
[`glm.cmp`](https://thomas-fung.github.io/mpcmp/reference/glm.cmp.md)
