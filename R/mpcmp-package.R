#' Mean-parametrized Conway-Maxwell-Poisson Regression
#'
#' @name mpcmp-package
#' @aliases mpcmp
#' @docType package
#' @title Mean-parametrized Conway-Maxwell-Poisson Regression
#' @keywords package
#' @references 
#' #' Fung, T., Alwan, A., Wishart, J. and Huang, A. (2019). The \code{mpcmp} package for 
#' Mean-parametrized Conway-Maxwell-Poisson Regression. 
#' 
#' Huang, A. (2017). Mean-parametrized Conway–Maxwell–Poisson regression models for 
#' dispersed counts. \emph{Statistical Modelling} \bold{17}, 359--380.
NULL

#' Attendance data set
#'
#' This data set gives the number of days absent 
#' from high school and the gender, maths score (standardized score out of 100) and 
#' academic programme (‘General’, ‘Academic’ and ‘Vocational’) of 314 students 
#' sampled from two urban high schools. The attendance data frame has 314 observations 
#' on 5 variables.
#'
#' @name attendance
#' @format A data frame with 314 observations on 5 variables.
#' \describe{
#' \item{id}{Identifier}
#' \item{gender}{gender}
#' \item{math}{standardized math score out of 100}
#' \item{daysabs}{number of days absent from high school}
#' \item{prog}{academic programme (‘General’, ‘Academic’ and ‘Vocational’)}
#' }
#' 
#' @docType data
#' @keywords datasets
#' @usage 
#' data(attendance)
#' @source \url{http://www.ats.ucla.edu/stat/stata/dae/nb_data.dta}
#' @examples 
#' ## For examples see example(glm.cmp)
NULL

#' Takeover Bids data set
#'
#' This data set gives the number of bids received by 126 US firms that were successful
#' targets of tender offers during the period 1978--1985, along with some explanatory 
#' variables on the defensive actions taken by management of target firm, firm-specific
#' characteristics and intervention taken by federal regulators. The \code{takeoverbids} 
#' data frame has 126 observations on 14 variables. The descriptions below are taken from
#'  Sáez-Castillo and Conde-Sánchez (2013).
#' 
#' 
#' @name takeoverbids
#' @format A data frame with 126 observations on 14 variables.
#' \describe{
#' \item{bidprem}{bid price divided by price 14 working days before bid}
#' \item{docno}{doc no}
#' \item{finrest}{indicator variable for proposed change in ownership structure}
#' \item{insthold}{percentage of stock held by institutions}
#' \item{leglrest}{indicator variable for legal defence by lawsuit}
#' \item{numbids}{number of bids recevied after the initial bid}
#' \item{obs}{Identifier}
#' \item{rearest}{indicator variable for proposed changes in asset structure}
#' \item{regulatn}{indicator variable for Department of Justice intervention}
#' \item{size}{total book value of assets in billions of dollars}
#' \item{takeover}{Indicator. 1 if the company was being taken over}
#' \item{weeks}{time in weeks between the initial and final offers}
#' \item{whtknght}{indicator varible for management invitation 
#' for friendly third-party bid}
#' \item{sizesq}{book value squared}
#' }
#' @docType data
#' @keywords datasets
#' @usage 
#' data(takeoverbids)
#' @references 
#' Cameron, A.C. and Johansson, P. (1997). Count Data Regression Models using Series 
#' Expansions: with Applications. \emph{Journal of Applied Econometrics} \bold{12} 203--223.
#' 
#' Cameron, A.C. and Trivedi P.K. (1998). Regression analysis of count data, Cambridge University Press, \url{http://cameron.econ.ucdavis.edu/racd/racddata.html} chapter 5.
#' 
#' Croissant Y (2011) Ecdat: Datasets for econometrics, R Package, version 0.1–6.1.
#' 
#' Jaggia, S. and Thosar, S. (1993). Multiple Bids as a Consequence of Target Management
#' Resistance \emph{Review of Quantitative Finance and Accounting} \bold{3}, 447--457.
#' 
#' @source
#' Journal of Applied Econometrics data archive: \url{http://qed.econ.queensu.ca/jae/}.
#' @examples 
#' ## For examples see example(glm.cmp)
NULL

#' Cotton Bolls data set
#'
#' This data set gives the observed number of bolls produced by the cotton 
#' plants at five growth stages: vegetative, flower-bud, blossom, fig and cotton boll; 
#' to examine the effect of five defoliation levels (0%, 25%, 50%, 75% and 100%). 
#' 
#' 
#' @name cottonbolls
#' @format A data frame with 125 observations on 4 variables.
#' \describe{
#' \item{nc}{number of bolls produced by two cotton plants at harvest}
#' \item{stages}{growth stage}
#' \item{def}{artificial defoliation level}
#' \item{def2}{square of def}
#' }
#' @docType data
#' @keywords datasets
#' @usage 
#' data(cottonbolls)
#' @references 
#' Zeviani, W.M., Riberio P.J. Jr., Bonat, W.H., Shimakura S.E. and Muniz J.A. (2014). 
#' The Gamma-count distribution in the analysis of experimental underdispersed data. 
#' \emph{Journal of Applied Statistics} \bold{41}, 2616--26.
#' 
#' @source
#' Supplementary Content of Zevini et al. (2014): 
#' \url{http://www.leg.ufpr.br/doku.php/publications:papercompanions:zeviani-jas2014}
#' 
#' @examples 
#' ## For examples see example(glm.cmp)
NULL

