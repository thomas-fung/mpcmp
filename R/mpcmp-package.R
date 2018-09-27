#' Mean-parametrized Conway-Maxwell Poisson Regression
#'
#' @name mpcmp-package
#' @docType package
#' @title Mean-parametrized Conway-Maxwell Poisson Regression
#' @keywords package
#' @references 
#' #' Fung, T., Alwan, A., Wishart, J. and Huang, A. (2018). The \code{mpcmp} package for 
#' Mean-parametrized Conway-Maxwell Poisson Regression. 
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
#' M.attendance = glm.cmp(daysabs ~ gender+math+prog, data=attendance)
#' @source \url{http://www.ats.ucla.edu/stat/stata/dae/nb_data.dta}
NULL

#' Takeover Bids data set
#'
#' This data set gives gives the number of bids received by 126 US firms that were successful
#' targets of tender offers during the period 1978--1985, along with some explanatory 
#' variables on the defensive actions taken by management of target firm, firm-specific
#' characteristics and any intervention taken by federal regulators. The \code{takeoverbids} 
#' data frame has 126 observations on 14 variables. 
#' 
#' 
#' @name takeoverbids
#' @format A data frame with 126 observations on 14 variables.
#' \describe{
#' \item{bidprem}{bid price divided by price 14 working days before bid}
#' \item{docno}{doc no.}
#' \item{finrest}{indicator variable for proposed change in ownership structure}
#' \item{insthold}{percentage of stock held by institutions}
#' \item{leglrest}{indicator variable for legal defence by lawsuit}
#' \item{numbids}{number of bids received}
#' \item{obs}{Identifier}
#' \item{rearest}{indicator variable for proposed changes in asset structure}
#' \item{regulartn}{indicator variable for Department of Justice intervention}
#' \item{size}{total book value of assets in billions of dollars}
#' \item{takeover}{Indicator. 1 if the company was being taken over}
#' \item{weeks}{}
#' \item{whtknght}{indicator varible for management invitation 
#' for friendly third-party bid}
#' }
#' @details Macropods defaecate randomly as they forage and scat 
#'   (faecal pellet) surveys are a reliable method for detecting the
#'   presence of rock-wallabies and other macropods. 
#'   Scats are used as an indication of spatial foraging patterns 
#'   of rock-wallabies and sympatric macropods. Scats deposited while
#'   foraging were not confused with scats deposited while
#'   resting because the daytime refuge areas of rock-wallabies
#'   were known in detail for each colony and no samples were
#'   taken from those areas. Each of the 200 sites were 
#'   examined separately to
#'   account for the different levels of predation risk and the
#'   abundance of rock-wallabies.
#' @docType data
#' @keywords datasets
#' @usage data(wallabies)
#' @references 
#'    Tuft KD, Crowther MS, Connell K, Mueller S and McArthur C (2011), 
#'    Predation risk and competitive interactions affect foraging of 
#'    an endangered refuge-dependent herbivore. Animal Conservation, 
#'    14: 447-457. doi: 10.1111/j.1469-1795.2011.00446.x
#' @examples
#' data(wallabies)
#' wdat = data.frame(subset(wallabies,select=-c(lat,long)), 
#'   EaD = wallabies$edible*wallabies$distance,
#'   EaS = wallabies$edible*wallabies$shelter,
#'   DaS = wallabies$distance*wallabies$shelter)
#' M1 = glm(rw~., family = binomial(link = "logit"), data = wdat)
NULL

