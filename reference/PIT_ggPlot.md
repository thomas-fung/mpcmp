# ggplot version of PIT Plots for a CMP Object

Two plots for the non-randomized PIT are currently available for
checking the distributional assumption of the fitted CMP model: the PIT
histogram, and the uniform Q-Q plot for PIT.

## Usage

``` r
gg_histcompPIT(
  object,
  bins = 10,
  ref_line = TRUE,
  col_line = "red",
  col_hist = "royal blue",
  size = 1
)

gg_qqcompPIT(
  object,
  bins = 10,
  col1 = "red",
  col2 = "#999999",
  lty1 = 1,
  lty2 = 2
)
```

## Arguments

- object:

  an object class "cmp", obtained from a call to `glm.cmp`.

- bins:

  numeric; the number of bins shown in the PIT histogram or the PIT Q-Q
  plot.

- ref_line:

  logical; if `TRUE` (default), the line for displaying the standard
  uniform distribution will be shown for the purpose of comparison.

- col_line:

  numeric or character: the colour of the reference line for comparison
  in PIT histogram.

- col_hist:

  numeric or character; the colour of the histogram for PIT.

- size:

  numeric; the line widths for the comparison line in PIT histogram.

- col1:

  numeric or character; the colour of the sample uniform Q-Q plot in
  PIT.

- col2:

  numeric or character; the colour of the theoretical uniform Q-Q plot
  in PIT.

- lty1:

  integer or character string: the line types for the sample uniform Q-Q
  plot in PIT, see ggplot2::linetype.

- lty2:

  an integer or character string: the line types for the theoretical
  uniform Q-Q plot in PIT, see ggplot2::linetype.

## Value

A `ggplot` object.

## Details

`histcompPIT` and `qqcompPIT`

The histogram and the Q-Q plot are used to compare the fitted profile
with a standard uniform distribution. If they match relatively well, it
means the CMP distribution is appropriate for the data.

The `histcompPIT` and `qqcompPIT` functions would provide the same two
plots but in base R format.

## References

Czado, C., Gneiting, T. and Held, L. (2009). Predictive model assessment
for count data. *Biometrics*, **65**, 1254–1261.

Dunsmuir, W.T.M. and Scott, D.J. (2015). The `glarma` Package for
Observation-Driven Time Series Regression of Counts. *Journal of
Statistical Software*, **67**, 1–36.

## See also

[`histcompPIT`](https://thomas-fung.github.io/mpcmp/reference/PIT_Plot.md),
[`qqcompPIT`](https://thomas-fung.github.io/mpcmp/reference/PIT_Plot.md),
[`plot.cmp`](https://thomas-fung.github.io/mpcmp/reference/plot.cmp.md)
and [`autoplot`](https://ggplot2.tidyverse.org/reference/autoplot.html).

## Examples

``` r
## For examples see example(autoplot)
```
