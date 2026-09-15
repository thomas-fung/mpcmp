# Update and Re-fit a COM-Poisson Model

`update` (i.e., `update.cmp`) will update and (by-default) re-fit a
model. It is identical to `update` in the `stats` package.

## Usage

``` r
# S3 method for class 'cmp'
update(object, formula., formula_nu., ..., evaluate = TRUE)
```

## Arguments

- object:

  an object class 'cmp', obtained from a call to `glm.cmp`.

- formula.:

  changes to the existing formula in `object` – see `update.formula`

- formula_nu.:

  changes to the existing formula_nu in `object` – see `update.formula`
  for details. It also accepts NULL to not regressing on the dispersion.

- ...:

  other arguments passed to or from other methods (currently unused).

- evaluate:

  logical; if `TRUE` evaluate the new call otherwise simply return the
  call

## Value

If `evaluate = TRUE` (the default), the re-fitted model object of class
`"cmp"`; otherwise the unevaluated call.

## See also

[`glm.cmp`](https://thomas-fung.github.io/mpcmp/reference/glm.cmp.md),
[`update.formula`](https://rdrr.io/r/stats/update.formula.html),
[`cmplrtest`](https://thomas-fung.github.io/mpcmp/reference/cmplrtest.md).

## Examples

``` r

# To update the mean regression formula
data(takeoverbids)

## Fit full model
M.bids.full <- glm.cmp(numbids ~ leglrest + rearest + finrest + whtknght
  + bidprem + insthold + size + sizesq + regulatn, data = takeoverbids)
M.bids.full
#> 
#> Call: glm.cmp(formula = numbids ~ leglrest + rearest + finrest + whtknght + 
#>     bidprem + insthold + size + sizesq + regulatn, data = takeoverbids)
#> 
#> Linear Model Coefficients:
#> (Intercept)     leglrest      rearest      finrest     whtknght      bidprem  
#>   0.9896300    0.2678800   -0.1731800    0.0677440    0.4812800   -0.6848200  
#>    insthold         size       sizesq     regulatn  
#>  -0.3678900    0.1793300   -0.0075823   -0.0375690  
#> 
#> Dispersion (nu): 1.75
#> Degrees of Freedom: 125 Total (i.e. Null);  116 Residual
#> Null Deviance: 182.3906 
#> Residual Deviance: 
#> AIC: 382.1753 
#> 

## Dropping whtknght
M.bids.null <- update(M.bids.full, . ~ . - whtknght)
M.bids.null
#> 
#> Call: glm.cmp(formula = numbids ~ leglrest + rearest + finrest + bidprem + 
#>     insthold + size + sizesq + regulatn, data = takeoverbids)
#> 
#> Linear Model Coefficients:
#> (Intercept)     leglrest      rearest      finrest      bidprem     insthold  
#>   1.3344000    0.3912400   -0.1464000    0.1342300   -0.7477300   -0.4406100  
#>        size       sizesq     regulatn  
#>   0.1878000   -0.0079426   -0.0445510  
#> 
#> Dispersion (nu): 1.54
#> Degrees of Freedom: 125 Total (i.e. Null);  117 Residual
#> Null Deviance: 164.7208 
#> Residual Deviance: 
#> AIC: 393.8799 
#> 

## To update the dispersion regression formula
data(sitophilus)

## Fit full model
M.sit.full <- glm.cmp(formula = ninsect ~ extract, formula_nu = ~extract, data = sitophilus)
M.sit.full
#> 
#> Call: glm.cmp(formula = ninsect ~ extract, formula_nu = ~extract, data = sitophilus)
#> 
#> Mean Model Coefficients:
#>   (Intercept)    extractLeaf  extractBranch    extractSeed  
#>       3.45000       -0.00637       -0.05210       -3.35000  
#> 
#> Dispersion Model Coefficients:
#>   (Intercept)    extractLeaf  extractBranch    extractSeed  
#>        -0.665         -0.383         -0.372         -0.118  
#> 
#> Degrees of Freedom: 39 Total (i.e. Null);  32 Residual
#> Null Deviance: 257.2108 
#> Residual Deviance: 
#> AIC: 260.8279 
#> 

## Dropping extract from the dispersion regression
M.sit.null1 <- update(M.sit.full, formula_nu = ~ . - extract)
M.sit.null1
#> 
#> Call: glm.cmp(formula = ninsect ~ extract, formula_nu = ~1, data = sitophilus)
#> 
#> Mean Model Coefficients:
#>   (Intercept)    extractLeaf  extractBranch    extractSeed  
#>       3.45000       -0.00637       -0.05210       -3.35000  
#> 
#> Dispersion Model Coefficients:
#> (Intercept)  
#>      -0.927  
#> 
#> Degrees of Freedom: 39 Total (i.e. Null);  35 Residual
#> Null Deviance: 234.4693 
#> Residual Deviance: 
#> AIC: 255.2668 
#> 

## To not regress on the dispersion at all
M.sit.null2 <- update(M.sit.full, formula_nu = NULL)
M.sit.null2
#> 
#> Call: glm.cmp(formula = ninsect ~ extract, data = sitophilus)
#> 
#> Linear Model Coefficients:
#>   (Intercept)    extractLeaf  extractBranch    extractSeed  
#>     3.4500000     -0.0063694     -0.0521290     -3.3547000  
#> 
#> Dispersion (nu): 0.396
#> Degrees of Freedom: 39 Total (i.e. Null);  36 Residual
#> Null Deviance: 234.4696 
#> Residual Deviance: 
#> AIC: 253.2668 
#> 
```
