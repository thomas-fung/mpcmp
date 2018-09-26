# mpcmp
Mean-parametrized Conway-Maxwell Poisson (mpcmp) Regression

The mpcmp package provides a collection of functions for estimation, testing and diagnostic checking for the mean-parametrized Conway-Maxwell Poisson regression model. 

We implement mean-parametrized COM-Poisson regreesion model (Huang (2017)) for under- and over-dispersed count data. The COM-Poisson currently only supports log-lienar mean models, however work is progressing to incorporate regression being linked to the dispersion parameter and a zero-inflated COM-Poisson model. 

You can see an example of the output here.

Installation

Check that you're running the most recent versions of your currently installed R packages:

update.packages()
Stable release on CRAN

The COMP package has been on CRAN since xxx 2018. You can install it from CRAN in the usual way:

install.packages("mpcmp")
library("mpcmp")

Development version on Github

You can use the devtools package to install the development version of mplot from GitHub:

# install.packages("devtools")
devtools::install_github("thomas-fung/mpcmp")
library(mpcmp)

Usage

A reference manual is available at thomas-htf.github.io/mpcmp

Citation

If you use this package , please use the following citation:

Fung, T, Alwan, A, Wishart, J and Huang, A. (2018). 

From R you can use:

citation("mpcmp")
toBibtex(citation("mpcmp"))
