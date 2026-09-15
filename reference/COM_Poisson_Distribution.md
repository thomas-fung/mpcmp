# The Conway-Maxwell-Poisson (COM-Poisson) Distribution.

Density, distribution function, quantile function and random generation
for the Conway-Maxwell-Poisson distribution with parameter `mu` and `nu`

## Usage

``` r
dcomp(
  x,
  mu,
  nu = 1,
  lambda,
  log.p = FALSE,
  lambdalb = 1e-10,
  lambdaub = 1000,
  maxlambdaiter = 1000,
  tol = 1e-06,
  summax
)

pcomp(
  q,
  mu,
  nu = 1,
  lambda,
  lower.tail = TRUE,
  log.p = FALSE,
  lambdalb = 1e-10,
  lambdaub = 1000,
  maxlambdaiter = 1000,
  tol = 1e-06,
  summax
)

qcomp(
  p,
  mu,
  nu = 1,
  lambda,
  lower.tail = TRUE,
  log.p = FALSE,
  lambdalb = 1e-10,
  lambdaub = 1000,
  maxlambdaiter = 1000,
  tol = 1e-06,
  summax
)

rcomp(
  n,
  mu,
  nu = 1,
  lambda,
  lambdalb = 1e-10,
  lambdaub = 1000,
  maxlambdaiter = 1000,
  tol = 1e-06,
  summax
)
```

## Arguments

- x, q:

  vector of quantiles

- mu, nu:

  mean and dispersion parameters. Must be strictly positive.

- lambda:

  an alternative way than mu to parametrized the distribution. Must be
  strictly positive

- log.p:

  logical; if `TRUE`, probabilities/densities \\p\\ are returned as
  \\log(p)\\.

- lambdalb, lambdaub:

  numeric: the lower and upper end points for the interval to be
  searched for lambda(s).

- maxlambdaiter:

  numeric: the maximum number of iterations allowed to solve for
  lambda(s).

- tol:

  numeric: the convergence threshold. A lambda is said to satisfy the
  mean constraint if the absolute difference between the calculated mean
  and mu is less than tol.

- summax:

  numeric; maximum number of terms to be considered in the truncated
  sum.

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\P(X \le x)\\,
  otherwise, \\P(X\>x)\\.

- p:

  vector of probabilities

- n:

  number of observations. If `length(n)` \> 1, the length is taken to be
  the number required.

## Value

`dcomp` gives the density, `pcomp` gives the distribution function,
`qcomp` gives the quantile function, and `rcomp` generates random
deviates.

Invalid arguments will result in return value `NaN`, with a warning.

The length of the results is determined by `n` for `rcomp`, and is the
maximum of the lengths of the numerical arguments for the other
functions.

The numerical arguments other than `n` are recycled to the length of the
results. Only the first argument of the logical arguments are used.

## Examples

``` r
dcomp(0:5, mu = 2, nu = 1.2)
#> [1] 0.11374731 0.27609781 0.29170831 0.18946331 0.08713145 0.03065722
pcomp(5, mu = 2, nu = 1.2)
#> [1] 0.9888054
p <- (1:9) / 10
qcomp(p, mu = 2, nu = 0.8)
#> [1] 0 1 1 1 2 2 3 3 4
rcomp(10, mu = 2, nu = 0.7)
#>  [1] 0 3 2 0 0 2 2 1 3 3
```
