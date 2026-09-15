# Likelihood Ratio Test for nu = 1 of a COM-Poisson model

Perform a likelihood ratio chi-squared test for nu = 1 of a COM-Poisson
model. The test statistics is calculated as *2\*(llik- llik_0)* where
*llik* and *llik_0* are the log-likelihood of a COM-Poisson and Poisson
model respectively. The test statistic has 1 degrees of freedom.

## Usage

``` r
LRTnu(object, digits = 3)
```

## Arguments

- object:

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

  the degrees of freedom for the test statistic (always 1).

- p.value:

  the p-value for the test.

- estimate:

  the log-likelihood of the fitted mean-CMP model and of the
  corresponding Poisson model.

- method:

  a character string describing the test.

- data.name:

  a character string giving the name of the model object.

Printing the returned object (e.g. via automatic printing at the
console) displays a formatted summary of the test.

## References

Huang, A. (2017). Mean-parametrized Conway-Maxwell-Poisson regression
models for dispersed counts. *Statistical Modelling* **17**, 359–380.

## Examples

``` r
data(takeoverbids)
M.bids <- glm.cmp(numbids ~ leglrest + rearest + finrest + whtknght
  + bidprem + insthold + size + sizesq + regulatn, data = takeoverbids)
LRTnu(M.bids)
#> 
#>  Likelihood ratio test for testing nu = 1
#> 
#> data:  M.bids
#> LR statistic = 9.72, df = 1, p-value = 0.001821
#> sample estimates:
#> log-lik for Mean-CMP  log-lik for Poisson 
#>                 -180                 -185 
#> 
```
