# Parameter Generator for nu

A function that use the arguments of a `glm.cmp` call to generate a
better initial `nu` estimate.

## Usage

``` r
getnu(
  param,
  y,
  xx,
  offset,
  llstart,
  fsscale = 1,
  lambdalb = 1e-10,
  lambdaub = 1000,
  maxlambdaiter = 1000,
  tol = 1e-06,
  summax = 100
)
```

## Arguments

- param:

  numeric vector: the model coefficients & the current value of `nu`. It
  is assumed that `nu` is in the last position of `param`.

- y:

  numeric vector: response variable

- xx:

  numeric matrix: the explanatory variables

- offset:

  numeric vector: a vector of length equal to the number of cases

- llstart:

  numeric: current log-likelihood value

- fsscale:

  numeric: a scaling factor (generally \>1) for the relaxed fisher
  scoring algorithm

- lambdalb, lambdaub:

  numeric: the lower and upper end points for the interval to be
  searched for lambda(s).

- maxlambdaiter:

  numeric: the maximum number of iterations allowed to solve for
  lambda(s).

- tol:

  numeric: the convergence threshold. A lambda is said to satisfy the
  mean constraint if the absolute difference between the calculated mean
  and a fitted values is less than tol.

- summax:

  maximum number of terms to be considered in the truncated sum

## Value

List containing the following:

- param:

  the model coefficients & the updated `nu`

- maxl:

  the updated log-likelihood

- fsscale:

  the final scaling factor used

## Details

From version 0.3.4, this function is no longer being used as part of the
estimation algorithm and this function will be defunct in our next
update.
