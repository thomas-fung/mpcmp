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
#>  [13]  0.09700635 -0.69146480  0.77969682  0.46652584  0.29965062  0.20915127
#>  [19]  0.65536346  0.36871504  1.78870981 -0.72277120 -0.22116010  0.68785292
#>  [25] -0.73335901 -0.54086034 -0.19484143  0.61935824  1.40171326 -0.34556604
#>  [31]  0.57976331 -0.61930564  1.66764210 -0.32674508 -1.19631776  3.00212218
#>  [37]  0.33943057  0.29955611 -0.05433707 -0.43315306 -0.03964149  0.21412266
#>  [43] -1.42418184  0.07067957  0.19904602 -0.42326501 -1.36287903 -0.82081648
#>  [49] -0.35647646  0.58104867 -0.02207807  0.72512952  0.40066704  0.21948704
#>  [55]  1.24537087 -0.94715388  0.57007228  0.19820112  0.35897952  0.55465718
#>  [61] -1.12223788  0.95901284 -2.57594638 -0.97647199  2.97257135 -0.21429119
#>  [67]  0.08282093 -0.24199884  0.69065204  0.09427416 -0.90708802  1.26858676
#>  [73]  1.06582648  2.19882669 -1.32972998  0.58666150  0.01491298 -0.32754280
#>  [79]  1.17220072 -0.17959482  0.53860595 -1.42781571 -0.75597401  0.32313015
#>  [85]  0.31710379  0.73548362  1.87904857  1.09434053 -1.27563878  0.08809889
#>  [91] -0.11024503 -0.33492048 -0.33475793 -1.10527210 -1.97614345 -0.61194972
#>  [97] -1.26431324 -0.29926741 -0.85858635  0.26496835  0.16016303 -1.04134637
#> [103] -0.66030831 -0.12961935  0.51059512  0.93994359 -1.48781668 -0.02032549
#> [109]  0.02041014 -0.82444946 -2.02222299 -0.05411675 -0.28020772  0.37053970
#> [115]  0.14630791 -0.90069021  0.55046992  2.90303429  0.09394954  1.13006033
#> [121]  0.78139897  1.93281653  0.88283234  1.07306945 -0.89773472 -2.10185581
#> 
#> $rtMid
#>   [1] -0.4898606052 -1.2821938866 -0.9046825631 -0.1188628661 -0.1341943952
#>   [6]  0.7653209654  0.2351273799  0.0171245480  0.3047524697 -1.1693426928
#>  [11]  0.0792278585  0.1389754525 -0.1034624751 -0.3370842293  0.3110988254
#>  [16]  0.1809578921  0.6063865740  0.1827372118  0.4748603628  0.1493556285
#>  [21]  2.0323888739 -0.4065244786 -0.0484976122  0.5966832639 -0.7121949353
#>  [26] -0.4681955935 -0.4067855841  0.0763118540  1.5778089677 -0.0213205174
#>  [31]  0.4153366142 -0.6181257783  1.5710150677 -0.4889168751 -0.6975868861
#>  [36]  2.8736704098  0.0293201149 -0.2299245554 -0.5074023913 -0.3107316170
#>  [41]  0.0040065546  0.3722261177 -1.0535891324  0.0433498533  0.2062825010
#>  [46] -0.6344033278 -1.0333900472 -1.0037282458  0.1400302408  0.6852283199
#>  [51] -0.3569787942  0.9064627311  0.3550496849  0.4009327343  1.1722251741
#>  [56] -0.5882093506  0.0362166245  0.1670724308  0.3196225946  0.7366028172
#>  [61] -1.0515786489  1.0488723319 -1.3465565763 -0.9436796455  3.0214268731
#>  [66]  0.1992771111  0.1892685799 -0.3062366880  0.3521605194  0.2239147663
#>  [71] -0.9333315833  1.2075756678  1.2234166316  2.2949753848 -1.3122964293
#>  [76]  0.1964248369  0.1269481758  0.0006044536  1.0115869879  0.1407696919
#>  [81] -0.0423125855 -1.0175105918 -0.7704861889 -0.0680272236  0.0658472526
#>  [86]  0.9119777110  1.8459551281  0.7570070437 -0.6610489089  0.1146635074
#>  [91] -0.4192417432 -0.6305074005 -0.2280802141 -0.7868264224 -1.8992452801
#>  [96] -0.7696441404 -1.2730186284 -0.4626315563 -0.8699211667 -0.1090825301
#> [101] -0.0729589351 -1.2752291746 -0.0714993128 -0.5068596334  0.0602518491
#> [106]  0.9397780473 -1.3557052945 -0.3086795425 -0.0736675539 -0.9880935055
#> [111] -2.2394275126 -0.1111248834  0.0010495725 -0.1014678556 -0.2458562634
#> [116] -0.6169651646  0.6335435909  2.7480914105  0.4121673593  1.3053444195
#> [121]  1.0357631163  1.8208308615  0.5919705287  0.6312364262 -0.9654756216
#> [126] -2.1524864287
#> 
```
