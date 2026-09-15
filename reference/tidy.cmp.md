# Tidy a(n) CMP model object

Tidy summarizes information about the components of a model. A model
component might be a single term in a regression, a single hypothesis, a
cluster, or a class. Exactly what tidy considers to be a model component
varies across models but is usually self-evident. If a model has several
distinct types of components, you will need to specify which components
to return.

## Usage

``` r
# S3 method for class 'cmp'
tidy(x, conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE, ...)
```

## Arguments

- x:

  an object class 'cmp' object, obtained from a call to `glm.cmp`

- conf.int:

  Logical indicating whether or not to include a confidence interval in
  the tidied output. Defaults to FALSE.

- conf.level:

  The confidence level to use for the confidence interval if conf.int =
  TRUE. Must be strictly greater than 0 and less than 1. Defaults to
  0.95, which corresponds to a 95 percent confidence interval.

- exponentiate:

  Logical indicating whether or not to exponentiate the the coefficient
  estimates.

- ...:

  other arguments passed to or from other methods (currently unused).

## Value

A
[`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
with columns:

- term:

  The name of the regression term.

- estimate:

  The estimated value of the regression term.

- std.error:

  The standard error of the regression term.

- statistic:

  The value of a test statistic to use in a hypothesis that the
  regression term is non-zero.

- p.value:

  The two-sided p-value associated with the observed statistic based on
  asymptotic normality.

- parameter:

  Only for varying dispersion models. Type of coefficient being
  estimated: 'mu', 'nu'

- conf.low:

  Lower bound on the confidence interval for the estimate.

- conf.high:

  Upper bound on the confidence interval for the estimate.

## Examples

``` r
data(attendance)
M.attendance <- glm.cmp(daysabs ~ gender + math + prog, data = attendance)
tidy(M.attendance)
#> # A tibble: 5 × 5
#>   term           estimate std.error statistic  p.value
#>   <chr>             <dbl>     <dbl>     <dbl>    <dbl>
#> 1 (Intercept)     2.71      0.190       14.3  4.05e-46
#> 2 gendermale     -0.215     0.117       -1.83 6.68e- 2
#> 3 math           -0.00632   0.00239     -2.65 8.04e- 3
#> 4 progAcademic   -0.425     0.170       -2.51 1.21e- 2
#> 5 progVocational -1.25      0.189       -6.62 3.65e-11
```
