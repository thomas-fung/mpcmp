# Summarizing COM-Poisson Model Fit

`summary` method for class `cmp`.

## Usage

``` r
# S3 method for class 'cmp'
summary(object, ...)

# S3 method for class 'summary.cmp'
print(
  x,
  digits = max(3, getOption("digits") - 3),
  signif.stars = getOption("show.signif.stars"),
  ...
)
```

## Arguments

- object:

  an object class 'cmp', obtained from a call to `glm.cmp`.

- ...:

  other arguments passed to or from other methods (currently unused).

- x:

  a result of the *default* method of
  [`summary()`](https://rdrr.io/r/base/summary.html).

- digits:

  numeric; minimum number of significant digits to be used for most
  numbers.

- signif.stars:

  logical. If TRUE, ‘significance stars’ are printed for each
  coefficient.

## Value

`summary.cmp` returns an object of class "summary.cmp", a list
containing at least the following components:

- call:

  the component from object.

- family:

  the component from object.

- deviance; residual_deviance:

  the component from object.

- df_residual:

  the component from object.

- df_null:

  the component from object.

- null_deviance:

  the component from object.

- deviance_resid:

  the deviance residuals: see residuals.cmp.

- coefficients:

  the matrix of coefficients, standard errors, z-values and p-values.

- df:

  a 3-vector of the rank of the model and the number of residual degrees
  of freedom, plus number of coefficients.

## Details

`print.summary.cmp` tries to be smart about formatting the coefficients,
standard errors and gives 'significance stars'. The `coefficients`
component of the result gives the estimated coefficients and their
estimated standard errors, together with their ratio. This third column
is labelled as `Z value` as the dispersion is fixed for this family. A
forth column gives the two-tailed p-value corresponding to `Z value`
based on the asymptotic Normal reference distribution.

## See also

[`coef.cmp`](https://thomas-fung.github.io/mpcmp/reference/coef.cmp.md),
[`fitted.cmp`](https://thomas-fung.github.io/mpcmp/reference/fitted.cmp.md),
[`glm.cmp`](https://thomas-fung.github.io/mpcmp/reference/glm.cmp.md).

## Examples

``` r
## For examples see example(glm.cmp)
```
