# Plot Diagnostic for a `glm.cmp` Object in ggplot style

`autoplot` uses ggplot2 to draw the diagnostic plots for a 'cmp' class
object. `gg_plot` is an *alias* for it.

## Usage

``` r
# S3 method for class 'cmp'
autoplot(
  object,
  which = c(1L, 2L, 6L, 8L),
  bins = 10,
  ask = TRUE,
  nrow = NULL,
  ncol = NULL,
  output_as_ggplot = TRUE,
  ...
)

gg_plot(
  object,
  which = c(1L, 2L, 6L, 8L),
  bins = 10,
  ask = TRUE,
  nrow = NULL,
  ncol = NULL,
  output_as_ggplot = TRUE,
  ...
)
```

## Arguments

- object:

  an object class 'cmp' object, obtained from a call to `glm.cmp`

- which:

  if a subset of plots is required, specify a subset of the numbers 1:8.
  See 'Details' below.

- bins:

  numeric; the number of bins shown in the PIT histogram or the PIT Q-Q
  plot.

- ask:

  logical; if `TRUE`, the user is asked before each plot.

- nrow:

  numeric; (optional) number of rows in the plot grid.

- ncol:

  numeric; (optional) number of columns in the plot grid.

- output_as_ggplot:

  logical; if `TRUE`, the function would return a list of `ggplot`
  objects; if `FALSE`, the function would return an `ggarrange` object.

- ...:

  other arguments passed to or from other methods (currently unused).

## Value

return a list of `ggplot` objects or a `ggarrange` object.

## Details

Eight plots (selectable by `which`) are currently available: a plot of
deviance residuals against fitted values, a non-randomized PIT
histogram, a uniform Q-Q plot for non-randomized PIT, a histogram of the
normal randomized residuals, a Q-Q plot of the normal randomized
residuals, a Scale-Location plot of sqrt(\| residuals \|) against fitted
values a plot of Cook's distances versus row labels a plot of pearson
residuals against leverage. By default, four plots (number 1, 2, 6, and
8 from this list of plots) are provided.

The 'Scale-Location' plot, also called 'Spread-Location' plot, takes the
square root of the absolute standardized deviance residuals
(*sqrt\|E\|*) in order to diminish skewness is much less skewed than
than *\|E\|* for Gaussian zero-mean E.

The 'Scale-Location' plot uses the standardized deviance residuals while
the Residual-Leverage plot uses the standardized pearson residuals. They
are given as \\R_i/\sqrt{1-h\_{ii}}\\ where \\h\_{ii}\\ are the diagonal
entries of the hat matrix.

The Residuals-Leverage plot shows contours of equal Cook's distance for
values of 0.5 and 1.

There are two plots based on the non-randomized probability integral
transformation (PIT) using
[`compPIT`](https://thomas-fung.github.io/mpcmp/reference/nrPIT.md).
These are a histogram and a uniform Q-Q plot. If the model assumption is
appropriate, these plots should reflect a sample obtained from a uniform
distribution.

There are also two plots based on the normal randomized residuals
calculated using
[`compnormRandPIT`](https://thomas-fung.github.io/mpcmp/reference/rPIT.md).
These are a histogram and a normal Q-Q plot. If the model assumption is
appropriate, these plots should reflect a sample obtained from a normal
distribution.

## See also

[`compPIT`](https://thomas-fung.github.io/mpcmp/reference/nrPIT.md),
[`compnormRandPIT`](https://thomas-fung.github.io/mpcmp/reference/rPIT.md),
[`glm.cmp`](https://thomas-fung.github.io/mpcmp/reference/glm.cmp.md)
and
[`plot.cmp`](https://thomas-fung.github.io/mpcmp/reference/plot.cmp.md).

## Examples

``` r
data(takeoverbids)
M.bids <- glm.cmp(numbids ~ leglrest + rearest + finrest + whtknght
  + bidprem + insthold + size + sizesq + regulatn, data = takeoverbids)

## The default plots are shown
gg_plot(M.bids) # or autoplot(M.bids)


## The plots for the non-randomized PIT
gg_plot(M.bids, which = c(2, 3)) # or autoplot(M.bids, which = c(2,3))
```
