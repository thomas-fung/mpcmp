# mpcmp: Mean-parametrized Conway-Maxwell Poisson Regression

# mpcmp <img src="man/figures/logo.svg" align ="right" alt="" width ="150"/>
<!-- badges: start -->
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/mpcmp)](https://cran.r-project.org/package=mpcmp)[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Travis build status](https://travis-ci.org/thomas-fung/mpcmp.svg?branch=master)](https://travis-ci.org/thomas-fung/mpcmp)
[![Codecov test coverage](https://codecov.io/gh/thomas-fung/mpcmp/branch/master/graph/badge.svg)](https://codecov.io/gh/thomas-fung/mpcmp?branch=master)
<!-- badges: end -->

The `mpcmp` package provides a collection of functions for estimation, testing and diagnostic checking for the mean-parametrized Conway-Maxwell Poisson (COM-Poisson) regression model for under- and over-dispersed count data of [Huang (2017)](https://doi.org/10.1177%2F1471082X17697749).

From version 0.3.0, `mpcmp` supports log-linear mean models, also allows one to incorporate regression being linked to the dispersion parameter.

Work is progressing to include a zero-inflated Conway-Maxwell-Poisson model. 

## Installation

### Stable release on CRAN

The ***mpcmp*** package has been on [CRAN](https://cran.r-project.org/package=mpcmp) since March 2019.  You can install it from CRAN in the usual way:

```s
install.packages("mpcmp")
library("mpcmp")
```

### Development version on Github

You can use the **devtools** package to install the development version of **mpcmp** from [GitHub](https://github.com/thomas-fung/mpcmp):

```s
# install.packages("devtools")
devtools::install_github("thomas-fung/mpcmp")
library(mpcmp)
```

## Usage

A reference manual is available at [thomas-fung.github.io/mpcmp](https://thomas-fung.github.io/mpcmp/)

## Citation

If you use this package to analyse your data, please use the following citation:

- Fung, T., Alwan, A., Wishart, J. and Huang, A. (2020). mpcmp: Mean-parametrized Conway-Maxwell Poisson Regression. R package version 0.3.0.

From R you can use:

```s
citation("mpcmp")
toBibtex(citation("mpcmp"))
```
