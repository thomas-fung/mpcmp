# Extract the Model Frame from a COM-Poisson Model Fit

An accessor function used to extract the model frame from a 'cmp'
object.

## Usage

``` r
# S3 method for class 'cmp'
model.frame(formula, ...)
```

## Arguments

- formula:

  an object class 'cmp' object, obtained from a call to `glm.cmp`

- ...:

  other arguments passed to or from other methods (currently unused).

## Value

The method will return the saved
[`data.frame`](https://rdrr.io/r/base/data.frame.html) used when fitting
the cmp model.

## See also

[`coef.cmp`](https://thomas-fung.github.io/mpcmp/reference/coef.cmp.md),
[`residuals.cmp`](https://thomas-fung.github.io/mpcmp/reference/residuals.cmp.md),
[`glm.cmp`](https://thomas-fung.github.io/mpcmp/reference/glm.cmp.md).
