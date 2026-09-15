# Extract Model Coefficients from a COM-Poisson Model Fit

An function used to extract model coefficients from a 'cmp' object.
`coefficients` is an alias for `coef`.

## Usage

``` r
# S3 method for class 'cmp'
coef(object, ...)
```

## Arguments

- object:

  an object class 'cmp' object, obtained from a call to `glm.cmp`

- ...:

  other arguments passed to or from other methods (currently unused).

## Value

Coefficients extracted from the object `object`.

## See also

[`fitted.cmp`](https://thomas-fung.github.io/mpcmp/reference/fitted.cmp.md),
[`residuals.cmp`](https://thomas-fung.github.io/mpcmp/reference/residuals.cmp.md),
[`glm.cmp`](https://thomas-fung.github.io/mpcmp/reference/glm.cmp.md).
