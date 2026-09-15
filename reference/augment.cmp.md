# Augment data with information from a(n) CMP model object

Augment accepts a model object and a dataset and adds information about
each observation in the dataset. Most commonly, this includes predicted
values in the .fitted column, residuals in the .resid column, and
standard errors for the fitted values in a .se.fit column. New columns
always begin with a . prefix to avoid overwriting columns in the
original dataset.

## Usage

``` r
# S3 method for class 'cmp'
augment(
  x,
  data = model.frame.cmp(x),
  newdata = NULL,
  type.predict = c("link", "response"),
  type.residuals = c("deviance", "pearson", "response"),
  se_fit = FALSE,
  ...
)
```

## Arguments

- x:

  an object class 'cmp' object, obtained from a call to
  [`glm.cmp`](https://thomas-fung.github.io/mpcmp/reference/glm.cmp.md)

- data:

  A base::data.frame or
  [`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
  containing the original data that was used to produce the object x.
  Defaults to model.frame.cmp(x) so that augment(my_fit) returns the
  augmented original data. **Do not** pass new data to the data
  argument. Augment will report information such as influence and cooks
  distance for data passed to the data argument. These measures are only
  defined for the original training data.

- newdata:

  A [`base::data.frame()`](https://rdrr.io/r/base/data.frame.html) or
  [`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
  containing all the original predictors used to create x. Defaults to
  NULL, indicating that nothing has been passed to newdata. If newdata
  is specified, the data argument will be ignored.

- type.predict:

  Passed to
  [`predict.cmp()`](https://thomas-fung.github.io/mpcmp/reference/predict.cmp.md)
  type argument. Defaults to `"link"`.

- type.residuals:

  Passed to
  [`residuals.cmp()`](https://thomas-fung.github.io/mpcmp/reference/residuals.cmp.md)
  type arguments. Defaults to `"deviance"`.

- se_fit:

  Logical indicating whether or not a .se.fit column should be added to
  the augmented output. Defaults to `FALSE`.

- ...:

  Additional arguments. Not used. Needed to match generic signature
  only. Cautionary note: Misspelled arguments will be absorbed in ...,
  where they will be ignored. If the misspelled argument has a default
  value, the default value will be used. For example, if you pass
  conf.level = 0.9, all computation will proceed using conf.level =
  0.95. Additionally, if you pass newdata = my_tibble to an augment()
  method that does not accept a newdata argument, it will use the
  default value for the data argument.

## Value

A
[`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
with columns:

- .cooksd:

  Cooks distance.

- .fitted:

  Fitted or predicted value.

- .hat:

  Diagonal of the hat matrix.

- .resid:

  The difference between observed and fitted values.

- .se.fit:

  Standard errors of fitted values.

- .sigma:

  Estimated residual standard deviation when corresponding observation
  is dropped from model.

- .std.resid:

  Standardised residuals.

## Details

Users may pass data to augment via either the data argument or the
newdata argument. If the user passes data to the data argument, it must
be exactly the data that was used to fit the model object. Pass datasets
to newdata to augment data that was not used during model fitting. This
still requires that all columns used to fit the model are present.

Augment will often behave differently depending on whether data or
newdata is given. This is because there is often information associated
with training observations (such as influences or related) measures that
is not meaningfully defined for new observations.

For convenience, many augment methods provide default data arguments, so
that augment(fit) will return the augmented training data. In these
cases, augment tries to reconstruct the original data based on the model
object with varying degrees of success.

The augmented dataset is always returned as a
[`tibble::tibble`](https://tibble.tidyverse.org/reference/tibble.html)
with the **same number of rows** as the passed dataset. This means that
the passed data must be coercible to a tibble.

We are in the process of defining behaviours for models fit with various
na.action arguments, but make no guarantees about behaviour when data is
missing at this time.

## Examples

``` r
data(attendance)
M.attendance <- glm.cmp(daysabs ~ gender + math + prog, data = attendance)
augment(M.attendance)
#> # A tibble: 314 × 9
#>    daysabs gender  math prog       .fitted .resid .std.resid    .hat   .cooksd
#>      <dbl> <fct>  <dbl> <fct>        <dbl>  <dbl>      <dbl>   <dbl>     <dbl>
#>  1       4 male      63 Academic     1.68  -0.264     -0.266 0.0115  0.000140 
#>  2       4 male      27 Academic     1.90  -0.464     -0.467 0.0106  0.000352 
#>  3       2 female    20 Academic     2.16  -1.14      -1.14  0.0128  0.00164  
#>  4       3 female    16 Academic     2.19  -0.907     -0.914 0.0140  0.00135  
#>  5       3 female     2 Academic     2.28  -0.976     -0.985 0.0201  0.00217  
#>  6      13 female    71 Academic     1.84   0.837      0.843 0.0137  0.00312  
#>  7      11 female    63 Academic     1.89   0.561      0.565 0.0117  0.00103  
#>  8       7 male       3 Academic     2.06  -0.109     -0.110 0.0184  0.0000424
#>  9      10 male      51 Academic     1.75   0.605      0.608 0.00968 0.00102  
#> 10       9 male      49 Vocational   0.936  1.47       1.48  0.0124  0.0124   
#> # ℹ 304 more rows
```
