# Solve for Lambda for a Particular Mean Parametrized COM-Poisson Distribution

Given a particular mean parametrized COM-Poisson distribution i.e. mu
and nu, this function is used to find a lambda that can satisfy the mean
constraint with a combination of bisection and Newton-Raphson updates.
The function is also vectorized but will only update those that have not
converged.

## Usage

``` r
comp_lambdas(
  mu,
  nu,
  lambdalb = 1e-10,
  lambdaub = 1000,
  maxlambdaiter = 1000,
  tol = 1e-06,
  lambdaint = 1,
  summax = 100
)

comp_lambdas_fixed_ub(
  mu,
  nu,
  lambdalb = 1e-10,
  lambdaub = 1000,
  maxlambdaiter = 1000,
  tol = 1e-06,
  lambdaint = 1,
  summax = 100
)
```

## Arguments

- mu, nu:

  mean and dispersion parameters. Must be straightly positive.

- lambdalb, lambdaub:

  numeric; the lower and upper end points for the interval to be
  searched for lambda(s).

- maxlambdaiter:

  numeric; the maximum number of iterations allowed to solve for
  lambda(s).

- tol:

  numeric; the convergence threshold. A lambda is said to satisfy the
  mean constraint if the absolute difference between the calculated mean
  and the corresponding mu values is less than tol.

- lambdaint:

  numeric vector; initial gauss for lambda(s).

- summax:

  maximum number of terms to be considered in the truncated sum

## Value

Both `comp_lambdas` and `comp_lambdas_fixed_ub` returns the lambda
value(s) that satisfies the mean constraint(s) as well as the current
lambdaub value. lambda value(s) returns by `comp_lambdas_fixed_ub` is
bounded by the lambdaub value. `comp_lambdas` has the extra ability to
scale up/down lambdaub to find the most appropriate lambda values.
