# Confidence Intervals for CMP Model Parameters

Computes confidence intervals for one or more parameters in a fitted
model.

## Usage

``` r
# S3 method for class 'cmp'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  an object class 'cmp', obtained from a call to
  [`glm.cmp`](https://thomas-fung.github.io/mpcmp/reference/glm.cmp.md).

- parm:

  a specification of which parameters are to be given confidence
  intervals, either a vector of numbers or a vector of names (comparing
  to those provided by [`coef()`](https://rdrr.io/r/stats/coef.html)) .
  If missing, all parameters are considered.

- level:

  the confidence level required.

- ...:

  other arguments passed to or from other methods (currently unused).

## Value

A matrix (or vector) with columns giving lower and upper confidence
limits for each parameter. These will be labelled as (1-level)/2 and 1 -
(1-level)/2 in % (by default 2.5% and 97.5%).

## Examples

``` r
data(attendance)
M.attendance <- glm.cmp(daysabs ~ gender + math + prog, data = attendance)
confint(M.attendance)
#>                       2.5%        97.5%
#> (Intercept)     2.34145406  3.087836298
#> gendermale     -0.44432495  0.014885329
#> math           -0.01099927 -0.001646776
#> progAcademic   -0.75758416 -0.093060290
#> progVocational -1.62526637 -0.882524629
confint(M.attendance, parm = "math", level = 0.9)
#>               5%          95%
#> math -0.01024745 -0.002398592
```
