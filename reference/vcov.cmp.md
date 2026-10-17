# Extracting the Variance-Covariance Matrix from a COM-Poisson Model Fit

Extracting the Variance-Covariance Matrix from a COM-Poisson Model Fit

## Usage

``` r
# S3 method for class 'cmp'
vcov(object, ...)
```

## Arguments

- object:

  an object class 'cmp' object, obtained from a call to `glm.cmp`

- ...:

  other arguments passed to or from other methods (currently unused).

## Value

The method will return the estimated covariances between the parameter
estimates of the fitted cmp model.

## Examples

``` r
data(attendance)
M.attendance <- glm.cmp(daysabs ~ gender + math + prog, data = attendance)
vcov(M.attendance)
#>                  (Intercept)    gendermale          math  progAcademic
#> (Intercept)     0.0362548810 -6.220523e-03 -2.500916e-04 -2.278627e-02
#> gendermale     -0.0062205234  1.372357e-02  5.319850e-06 -7.355542e-04
#> math           -0.0002500916  5.319850e-06  5.692439e-06  1.039728e-05
#> progAcademic   -0.0227862685 -7.355542e-04  1.039728e-05  2.873856e-02
#> progVocational -0.0189576491 -5.111475e-04 -7.981398e-05  2.253265e-02
#>                progVocational
#> (Intercept)     -1.895765e-02
#> gendermale      -5.111475e-04
#> math            -7.981398e-05
#> progAcademic     2.253265e-02
#> progVocational   3.590207e-02
```
