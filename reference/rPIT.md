# Random Normal Probability Integral Transform

A function to create the normal conditional (randomized) quantile
residuals. The majority of the code and descriptions are taken from
Dunsmuir and Scott (2015).

## Usage

``` r
compnormRandPIT(object)
```

## Arguments

- object:

  an object class "cmp", obtained from a call to `glm.cmp`.

## Value

A list consisting of two elements:

- rt:

  the normal conditional randomized quantile residuals

- rdMid:

  the midpoints of the predictive probability intervals

## Details

The function `compPredProb` produces the non-randomized probability
integral transform(PIT). It returns estimates of the cumulative
predictive probabilities as upper and lower bounds of a collection of
intervals. If the model is correct, a histogram drawn using these
estimated probabilities should resemble a histogram obtained from a
sample from the uniform distribution.

This function aims to produce observations which instead resemble a
sample from a normal distribution. Such a sample can then be examined by
the usual tools for checking normality, such as histograms and normal
Q-Q plots.

For each of the intervals produced by `compnormRandPIT`, a random
uniform observation is generated, which is then converted to a normal
observation by applying the inverse standard normal distribution
function (using `qnorm`). The vector of these values is returned by the
function in the list element `rt`. In addition non-random observations
which should appear similar to a sample from a normal distribution are
obtained by applying `qnorm` to the mid-points of the predictive
distribution intervals. The vector of these values is returned by the
function in the list element `rtMid`.

## References

Berkowitz, J. (2001). Testing density forecasts, with applications to
risk management. *Journal of Business and Economic Statistics*, **19**,
465–474.

Dunn, P. K. and Smyth, G. K. (1996). Randomized quantile residuals.
*Journal of Computational and Graphical Statistics*, **5**, 236–244.

Dunsmuir, W.T.M. and Scott, D.J. (2015). The `glarma` Package for
Observation-Driven Time Series Regression of Counts. *Journal of
Statistical Software*, **67**, 1–36.

## Examples

``` r
data(takeoverbids)
M.bids <- glm.cmp(numbids ~ leglrest + rearest + finrest + whtknght
  + bidprem + insthold + size + sizesq + regulatn, data = takeoverbids)
compnormRandPIT(M.bids)
#> $rt
#>   [1] -0.20016809 -1.81257497 -1.51743440 -0.32208532 -0.24345221  0.55642354
#>   [7]  0.15153776 -0.49321625  0.17991618 -0.71811993 -0.15837151  0.29653545
#>  [13]  0.09700635 -0.69146480  0.77969682  0.46652584  0.29965061  0.20915127
#>  [19]  0.65536346  0.36871504  1.78870981 -0.72277120 -0.22116010  0.68785292
#>  [25] -0.73335901 -0.54086034 -0.19484143  0.61935824  1.40171326 -0.34556604
#>  [31]  0.57976331 -0.61930564  1.66764210 -0.32674508 -1.19631776  3.00212218
#>  [37]  0.33943057  0.29955611 -0.05433707 -0.43315306 -0.03964149  0.21412266
#>  [43] -1.42418185  0.07067957  0.19904602 -0.42326501 -1.36287903 -0.82081648
#>  [49] -0.35647646  0.58104867 -0.02207807  0.72512952  0.40066704  0.21948704
#>  [55]  1.24537087 -0.94715388  0.57007228  0.19820112  0.35897952  0.55465718
#>  [61] -1.12223788  0.95901284 -2.57594638 -0.97647200  2.97257135 -0.21429119
#>  [67]  0.08282093 -0.24199884  0.69065204  0.09427416 -0.90708802  1.26858676
#>  [73]  1.06582648  2.19882669 -1.32972998  0.58666150  0.01491298 -0.32754280
#>  [79]  1.17220072 -0.17959482  0.53860595 -1.42781571 -0.75597401  0.32313015
#>  [85]  0.31710379  0.73548362  1.87904857  1.09434053 -1.27563878  0.08809889
#>  [91] -0.11024503 -0.33492048 -0.33475793 -1.10527210 -1.97614345 -0.61194972
#>  [97] -1.26431324 -0.29926741 -0.85858635  0.26496835  0.16016303 -1.04134637
#> [103] -0.66030831 -0.12961935  0.51059512  0.93994359 -1.48781668 -0.02032549
#> [109]  0.02041014 -0.82444946 -2.02222299 -0.05411675 -0.28020772  0.37053970
#> [115]  0.14630791 -0.90069021  0.55046992  2.90303429  0.09394953  1.13006033
#> [121]  0.78139897  1.93281653  0.88283234  1.07306945 -0.89773472 -2.10185581
#> 
#> $rtMid
#>   [1] -0.4898606057 -1.2821938867 -0.9046825634 -0.1188628663 -0.1341943954
#>   [6]  0.7653209649  0.2351273795  0.0171245478  0.3047524696 -1.1693426929
#>  [11]  0.0792278584  0.1389754521 -0.1034624756 -0.3370842295  0.3110988250
#>  [16]  0.1809578920  0.6063865735  0.1827372114  0.4748603624  0.1493556284
#>  [21]  2.0323888736 -0.4065244788 -0.0484976126  0.5966832630 -0.7121949355
#>  [26] -0.4681955938 -0.4067855844  0.0763118539  1.5778089671 -0.0213205176
#>  [31]  0.4153366138 -0.6181257786  1.5710150671 -0.4889168754 -0.6975868864
#>  [36]  2.8736704081  0.0293201148 -0.2299245556 -0.5074023916 -0.3107316172
#>  [41]  0.0040065545  0.3722261174 -1.0535891327  0.0433498529  0.2062825009
#>  [46] -0.6344033284 -1.0333900475 -1.0037282461  0.1400302406  0.6852283194
#>  [51] -0.3569787944  0.9064627309  0.3550496840  0.4009327342  1.1722251736
#>  [56] -0.5882093509  0.0362166244  0.1670724301  0.3196225942  0.7366028169
#>  [61] -1.0515786490  1.0488723314 -1.3465565764 -0.9436796458  3.0214268723
#>  [66]  0.1992771107  0.1892685795 -0.3062366882  0.3521605190  0.2239147662
#>  [71] -0.9333315837  1.2075756668  1.2234166311  2.2949753838 -1.3122964298
#>  [76]  0.1964248365  0.1269481756  0.0006044535  1.0115869876  0.1407696912
#>  [81] -0.0423125857 -1.0175105921 -0.7704861895 -0.0680272238  0.0658472525
#>  [86]  0.9119777104  1.8459551277  0.7570070434 -0.6610489092  0.1146635070
#>  [91] -0.4192417434 -0.6305074008 -0.2280802143 -0.7868264227 -1.8992452803
#>  [96] -0.7696441407 -1.2730186285 -0.4626315566 -0.8699211671 -0.1090825303
#> [101] -0.0729589353 -1.2752291754 -0.0714993130 -0.5068596336  0.0602518486
#> [106]  0.9397780470 -1.3557052949 -0.3086795427 -0.0736675540 -0.9880935056
#> [111] -2.2394275128 -0.1111248838  0.0010495723 -0.1014678557 -0.2458562636
#> [116] -0.6169651649  0.6335435906  2.7480914094  0.4121673590  1.3053444188
#> [121]  1.0357631161  1.8208308609  0.5919705281  0.6312364259 -0.9654756216
#> [126] -2.1524864294
#> 
```
