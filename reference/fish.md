# Fish data set

This data set gives the the number of fish species in lakes of the
world; to examine the effect of the surface area of the lakes. The
latitude of the lakes are also recorded.

## Usage

``` r
data(fish)
```

## Format

A data frame with 70 observations on 4 variables.

- lake:

  name of the lakes

- species:

  number of fish species in lakes

- area:

  surface area (km squared)

- latitude:

  latitude of the lakes

## Value

A data frame with 70 rows and 4 variables.

## Details

This data set is also used to illustrate that the fitting algorithm can
handle some larger count data.

## References

Barbour, C. D. and Brown, J. H. (1974). Fish species diversity in lakes.
*The American Naturalist*, **108**, 473–488.

## Examples

``` r
### Barbour & Brown (1974): Overdispersed Fish data
# \donttest{
data(fish)
M.fish <- glm.cmp(species ~ 1 + log(area), data = fish)
M.fish
#> 
#> Call: glm.cmp(formula = species ~ 1 + log(area), data = fish)
#> 
#> Linear Model Coefficients:
#> (Intercept)    log(area)  
#>     2.32810      0.18286  
#> 
#> Dispersion (nu): 0.0184
#> Degrees of Freedom: 69 Total (i.e. Null);  68 Residual
#> Null Deviance: 101.6675 
#> Residual Deviance: 59.45133 
#> AIC: 638.8532 
#> 
summary(M.fish)
#> 
#> Call: glm.cmp(formula = species ~ 1 + log(area), data = fish)
#> 
#> Deviance Residuals: 
#>     Min       1Q   Median       3Q      Max  
#> -2.1215  -0.8038  -0.3327   0.4276   2.6966  
#> 
#> Linear Model Coefficients:
#>             Estimate Std.Err Z value Pr(>|z|)    
#> (Intercept)  2.32810 0.23410   9.945  < 2e-16 ***
#> log(area)    0.18286 0.02872   6.368 1.91e-10 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for Mean-CMP estimated to be 0.01841)
#> 
#> 
#>     Null deviance: 101.668  on 69 degrees of freedom
#> Residual deviance:  59.451  on 68 degrees of freedom
#> 
#> AIC: 638.8532 
#> 
# }
```
