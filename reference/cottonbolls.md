# Cotton Bolls data set

This data set gives the observed number of bolls produced by the cotton
plants at five growth stages: vegetative, flower-bud, blossom, fig and
cotton boll; to examine the effect of five defoliation levels (0\\

## Usage

``` r
data(cottonbolls)
```

## Format

A data frame with 125 observations on 4 variables.

- nc:

  number of bolls produced by two cotton plants at harvest

- stages:

  growth stage

- def:

  artificial defoliation level

- def2:

  square of def

## Source

Supplementary Content of Zeviani et al. (2014):
<http://www.leg.ufpr.br/doku.php/publications:papercompanions:zeviani-jas2014>

## Value

A data frame with 125 rows and 4 variables.

## References

Zeviani, W.M., Riberio P.J. Jr., Bonat, W.H., Shimakura S.E. and Muniz
J.A. (2014). The Gamma-count distribution in the analysis of
experimental underdispersed data. *Journal of Applied Statistics*
**41**, 2616–26.

## Examples

``` r
### Huang (2017) Page 373--375: Underdispersed Cotton bolls data
### Model fitting for predictor V

data(cottonbolls)
M.bolls <- glm.cmp(nc ~ 1 + stages:def + stages:def2, data = cottonbolls)
M.bolls
#> 
#> Call: glm.cmp(formula = nc ~ 1 + stages:def + stages:def2, data = cottonbolls)
#> 
#> Linear Model Coefficients:
#>            (Intercept)    stagesvegetative:def    stagesflower bud:def  
#>              2.1896000               0.4368500               0.2897100  
#>      stagesblossom:def           stagesfig:def   stagescotton boll:def  
#>             -1.2425000               0.3648200               0.0089438  
#>  stagesvegetative:def2   stagesflower bud:def2      stagesblossom:def2  
#>             -0.8052100              -0.4878400               0.6728600  
#>         stagesfig:def2  stagescotton boll:def2  
#>             -1.3103000              -0.0199660  
#> 
#> Dispersion (nu): 4.85
#> Degrees of Freedom: 124 Total (i.e. Null);  114 Residual
#> Null Deviance: 345.9363 
#> Residual Deviance: 
#> AIC: 440.8229 
#> 
summary(M.bolls)
#> 
#> Call: glm.cmp(formula = nc ~ 1 + stages:def + stages:def2, data = cottonbolls)
#> 
#> Deviance Residuals: 
#>      Min        1Q    Median        3Q       Max  
#> -2.29766  -0.68296  -0.04972   0.74987   2.74460  
#> 
#> Linear Model Coefficients:
#>                         Estimate   Std.Err Z value Pr(>|z|)    
#> (Intercept)             2.189562  0.029407  74.456  < 2e-16 ***
#> stagesvegetative:def    0.436854  0.239663   1.823  0.06834 .  
#> stagesflower bud:def    0.289708  0.235822   1.229  0.21926    
#> stagesblossom:def      -1.242507  0.283166  -4.388 1.14e-05 ***
#> stagesfig:def           0.364816  0.264720   1.378  0.16817    
#> stagescotton boll:def   0.008944  0.233866   0.038  0.96949    
#> stagesvegetative:def2  -0.805212  0.271949  -2.961  0.00307 ** 
#> stagesflower bud:def2  -0.487844  0.263290  -1.853  0.06390 .  
#> stagesblossom:def2      0.672864  0.319551   2.106  0.03523 *  
#> stagesfig:def2         -1.310267  0.316982  -4.134 3.57e-05 ***
#> stagescotton boll:def2 -0.019966  0.256498  -0.078  0.93795    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for Mean-CMP estimated to be 4.852)
#> 
#> 
#>     Null deviance: 345.94  on 124 degrees of freedom
#> Residual deviance: 125.26  on 114 degrees of freedom
#> 
#> AIC: 440.8229 
#> 
```
