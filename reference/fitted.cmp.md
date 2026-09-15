# Extract Fitted Values from a COM-Poisson Model Fit

An accessor function used to extract the fitted values from a 'cmp'
object. `fitted.values` is an alias for `fitted`.

## Usage

``` r
# S3 method for class 'cmp'
fitted(object, ...)
```

## Arguments

- object:

  an object class 'cmp' object, obtained from a call to `glm.cmp`

- ...:

  other arguments passed to or from other methods (currently unused).

## Value

Fitted values `mu` extracted from the object `object`.

## See also

[`coef.cmp`](https://thomas-fung.github.io/mpcmp/reference/coef.cmp.md),
[`residuals.cmp`](https://thomas-fung.github.io/mpcmp/reference/residuals.cmp.md),
[`glm.cmp`](https://thomas-fung.github.io/mpcmp/reference/glm.cmp.md).
