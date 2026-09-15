# Print Method for a `cmp_autoplot_list` Object

Re-displays the combined diagnostic plot produced by
[`autoplot.cmp`](https://thomas-fung.github.io/mpcmp/reference/autoplot.cmp.md)
when `output_as_ggplot = TRUE`. This avoids the default list-print
behaviour (`[[1]]`, `[[2]]`, ... headers with each plot printed
separately) when the list of `ggplot` objects returned by
`autoplot.cmp`/`gg_plot` is printed again, e.g. after being assigned to
a variable.

## Usage

``` r
# S3 method for class 'cmp_autoplot_list'
print(x, ...)
```

## Arguments

- x:

  an object of class `cmp_autoplot_list`, as returned by
  [`autoplot.cmp`](https://thomas-fung.github.io/mpcmp/reference/autoplot.cmp.md)
  with `output_as_ggplot = TRUE`.

- ...:

  other arguments passed to or from other methods (currently unused).

## Value

`x` is returned invisibly; `print.cmp_autoplot_list` is called for its
side effect of drawing the combined diagnostic plot on the current
graphics device.
