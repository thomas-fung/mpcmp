# mpcmp: Mean-parametrized Conway-Maxwell Poisson Regression

The `mpcmp` package provides a collection of functions for estimation, testing and diagnostic checking for the mean-parametrized Conway-Maxwell Poisson (COM-Poisson) regression model for under- and over-dispersed count data of [Huang (2017)](https://doi.org/10.1177%2F1471082X17697749).

The `mpcmp` currently only supports log-lienar mean models, however work is progressing to incorporate regression being linked to the dispersion parameter and a zero-inflated Conway-Maxwell Poisson model. 

## Installation

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

- Fung, T., Alwan, A., Wishart, J. and Huang, A. (2019). "The mpcmp package for Mean-parametrized Conway-Maxwell Poisson Regression."

From R you can use:

```s
citation("mpcmp")
toBibtex(citation("mpcmp"))
```
