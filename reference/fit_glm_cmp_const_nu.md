# Fit a Mean Parametrized Conway-Maxwell Poisson Generalized Linear Model with constant dispersion.

This is a workhorse function in which glm.cmp to call upon to fit a
mean-parametrized Conway-Maxwell Poisson generalized linear model with
constant dispersion.

## Usage

``` r
fit_glm_cmp_const_nu(
  y = y,
  X = X,
  offset = offset,
  betastart = betastart,
  lambdalb = lambdalb,
  lambdaub = lambdaub,
  maxlambdaiter = maxlambdaiter,
  tol = tol
)
```

## Arguments

- y:

  the response y vector.

- X:

  the design matrix for regressing the mean

- offset, :

  this can be used to specify an *a priori* known component to be
  included in the linear predictor for mean during fitting. This should
  be `NULL` or a numeric vector

- betastart:

  starting values for the parameters in the linear predictor for mu.

- lambdalb, lambdaub:

  numeric: the lower and upper end points for the interval to be
  searched for lambda(s). The default value for lambdaub should be
  sufficient for small to moderate size nu. If nu is large and required
  a larger `lambdaub`, the algorithm will scale up `lambdaub`
  accordingly.

- maxlambdaiter:

  numeric: the maximum number of iterations allowed to solve for
  lambda(s).

- tol:

  numeric: the convergence threshold. A lambda is said to satisfy the
  mean constraint if the absolute difference between the calculated mean
  and a fitted values is less than tol.

## Value

A fitted model object of class `cmp` similar to one obtained from `glm`
or `glm.nb`.

## Examples

``` r
## For examples see example(glm.cmp)
```
