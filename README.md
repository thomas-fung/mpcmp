# mpcmp: Mean-parametrized Conway-Maxwell Poisson Regression <img src="man/figures/logo.svg" align="right" alt="" width="150"/>

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/mpcmp)](https://cran.r-project.org/package=mpcmp)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/thomas-fung/mpcmp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/thomas-fung/mpcmp/actions/workflows/R-CMD-check.yaml)
[![R-hub](https://github.com/thomas-fung/mpcmp/actions/workflows/rhub.yaml/badge.svg)](https://github.com/thomas-fung/mpcmp/actions/workflows/rhub.yaml)
[![Codecov test coverage](https://codecov.io/gh/thomas-fung/mpcmp/graph/badge.svg)](https://app.codecov.io/gh/thomas-fung/mpcmp)
[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/thomas-fung/mpcmp)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/mpcmp)](https://cran.r-project.org/package=mpcmp)
[![pkgdown](https://img.shields.io/badge/pkgdown-docs-blue.svg)](https://thomas-fung.github.io/mpcmp/)
<!-- badges: end -->

The `mpcmp` package provides a collection of functions for estimation, testing and diagnostic checking for the mean-parametrized Conway-Maxwell Poisson (COM-Poisson) regression model for under- and over-dispersed count data of [Huang (2017)](https://doi.org/10.1177/2F1471082X17697749).

From version 0.3.0, `mpcmp` supports log-linear mean models, also allows one to incorporate regression being linked to the dispersion parameter.

Work is progressing to include a zero-inflated Conway-Maxwell-Poisson model. 

## Installation

### Stable release on CRAN

The ***mpcmp*** package has been on [CRAN](https://cran.r-project.org/package=mpcmp) since March 2019.  You can install it from CRAN in the usual way:

```s
# install.packages("pak")
pak::pak("mpcmp")
```

### Development version on Github

You can use the **devtools** package to install the development version of **mpcmp** from [GitHub](https://github.com/thomas-fung/mpcmp):

```s
# install.packages("pak")
pak::pak("thomas-fung/mpcmp")
library(mpcmp)
```

## Usage

A reference manual and some examples are available at [thomas-fung.github.io/mpcmp](https://thomas-fung.github.io/mpcmp/)

## Citation

If you use this package to analyse your data, please use the following citation:

- Fung, T., Alwan, A., Wishart, J. and Huang, A. (2026). mpcmp: Mean-parametrized Conway-Maxwell Poisson Regression. R package version 0.4.0.

From R you can use:

```s
citation("mpcmp")
toBibtex(citation("mpcmp"))
```
