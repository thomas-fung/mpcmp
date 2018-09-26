# mpcmp: Mean-parametrized Conway-Maxwell Poisson Regression

The `mpcmp` package provides a collection of functions for estimation, testing and diagnostic checking for the mean-parametrized Conway-Maxwell Poisson (COM-Poisson) regression model for under- and over-dispersed count data of [Huang (2017)](https://doi.org/10.1177%2F1471082X17697749).

The `mpcmp` currently only supports log-lienar mean models, however work is progressing to incorporate regression being linked to the dispersion parameter and a zero-inflated Conway-Maxwell Poisson model. 

## Installation

### Stable release on CRAN

The mpcmp package has been on CRAN since xxx 2018. You can install it from CRAN in the usual way:

```s
install.packages("mpcmp")
library("mpcmp")
```

### Development version on Github

You can use the **devtools** package to install the development version of **mpcmp** from GitHub:

```s
# install.packages("devtools")
devtools::install_github("thomas-fung/mpcmp")
library(mpcmp)
```

## Usage

A reference manual is available at thomas-htf.github.io/mpcmp

## Citation

If you use this package to analyse your data, please use the following citation:

- Fung, T., Alwan, A., Wishart, J. and Huang, A. (2018). "The mpcmp package for Mean-parametrized Conway-Maxwell Poisson Regression."

From R you can use:

```s
citation("mpcmp")
toBibtex(citation("mpcmp"))
```