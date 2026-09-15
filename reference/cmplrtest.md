# Likelihood Ratio Test for nested COM-Poisson models

Perform a likelihood ratio chi-squared test between nested COM-Poisson
models. The test statistics is calculated as *2\*(llik- llik_0)*. The
test statistics has degrees of freedom *r* where *r* is the difference
in the number of parameters between the full and null models.

## Usage

``` r
cmplrtest(object1, object2, digits = 3)
```

## Arguments

- object1:

  an object class 'cmp', obtained from a call to `glm.cmp`

- object2:

  an object class 'cmp', obtained from a call to `glm.cmp`

- digits:

  numeric; minimum number of significant digits to be used for most
  numbers.

## Value

An object of class `"htest"` (see
[`t.test`](https://rdrr.io/r/stats/t.test.html) for the generic
structure), with components:

- statistic:

  the likelihood ratio test statistic.

- parameter:

  the degrees of freedom for the test statistic.

- p.value:

  the p-value for the test.

- method:

  a character string describing the test.

- data.name:

  a character string giving the names of the two model objects compared.

Printing the returned object (e.g. via automatic printing at the
console) displays a formatted summary of the test.

## References

Huang, A. (2017). Mean-parametrized Conway-Maxwell-Poisson regression
models for dispersed counts. *Statistical Modelling* **17**, 359–380.

## See also

[`glm.cmp`](https://thomas-fung.github.io/mpcmp/reference/glm.cmp.md),
[`update.cmp`](https://thomas-fung.github.io/mpcmp/reference/update.cmp.md)

## Examples

``` r

## Testing for the mean coefficients
data(takeoverbids)

## Fit full model
M.bids.full <- glm.cmp(numbids ~ leglrest + rearest + finrest + whtknght
  + bidprem + insthold + size + sizesq + regulatn, data = takeoverbids)

## Fit null model; without whtknght
M.bids.null <- update(M.bids.full, . ~ . - whtknght)

## Likelihood ratio test for the nested models
cmplrtest(M.bids.full, M.bids.null) # order of objects is not important
#> 
#>  Likelihood ratio test for testing both COM-Poisson models are
#>  equivalent
#> 
#> data:  M.bids.full vs. M.bids.null
#> LR statistic = 13.7, df = 1, p-value = 0.0002139
#> 

## Testing for dispersion coefficients
data(sitophilus)
M.sit.full <- glm.cmp(formula = ninsect ~ extract, formula_nu = ~extract, data = sitophilus)

## Fit null model; dropping extract from dispersion equation
M.sit.null1 <- update(M.sit.full, formula_nu. = ~1)
cmplrtest(M.sit.null1, M.sit.full)
#> 
#>  Likelihood ratio test for testing both COM-Poisson models are
#>  equivalent
#> 
#> data:  M.sit.null1 vs. M.sit.full
#> LR statistic = 0.439, df = 3, p-value = 0.9321
#> 

## Fit null model; using constant dispersion specification
M.sit.null2 <- update(M.sit.full, formula_nu. = NULL)
cmplrtest(M.sit.null2, M.sit.full)
#> 
#>  Likelihood ratio test for testing both COM-Poisson models are
#>  equivalent
#> 
#> data:  M.sit.null2 vs. M.sit.full
#> LR statistic = 0.439, df = 3, p-value = 0.9321
#> 
```
