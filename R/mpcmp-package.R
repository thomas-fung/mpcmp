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


#' Blood and other measurements in diabetics
#'
#' The diabetes data frame has 442 rows and 11 columns.
#' These are the data used in Efron et al. (2004).
#'
#' @name diabetes
#' @format A data frame with 442 observations on 11 variables.
#' \describe{
#' \item{age}{Age}
#' \item{sex}{Gender}
#' \item{bmi}{Body mass index}
#' \item{map}{Mean arterial pressure (average blood pressure)}
#' \item{tc}{Total cholesterol (mg/dL)? Desirable range: below 200 mg/dL}
#' \item{ldl}{Low-density lipoprotein ("bad" cholesterol)? 
#'            Desirable range: below 130 mg/dL }
#' \item{hdl}{High-density lipoprotein ("good" cholesterol)? 
#'            Desirable range: above 40 mg/dL}
#' \item{tch}{Blood serum measurement}
#' \item{ltg}{Blood serum measurement}
#' \item{glu}{Blood serum measurement (glucose?)}
#' \item{y}{A quantitative measure of disease progression 
#'          one year after baseline}
#' }
#' @details Data sourced from http://web.stanford.edu/~hastie/Papers/LARS
#' @docType data
#' @keywords datasets
#' @usage data(diabetes)
#' @references Efron, B., Hastie, T., Johnstone, I., Tibshirani, R., (2004).
#'   Least angle regression. The Annals of Statistics 32(2) 407-499.
#'   DOI: 10.1214/009053604000000067
#' @examples
#' data(diabetes)
#' full.mod = lm(y~.,data=diabetes)
NULL




#' Artificial example
#'
#' An artificial data set which causes stepwise regression
#' procedures to select a non-parsimonious model.
#' The true model is a simple linear regression of
#' y against x8.
#'
#' @name artificialeg
#' @format A data frame with 50 observations on 10 variables.
#' @details Inspired by the pathoeg data set in the MPV pacakge.
#' @docType data
#' @keywords datasets
#' @usage data(artificialeg)
#' @examples
#' data(artificialeg)
#' full.mod = lm(y~.,data=artificialeg)
#' step(full.mod)
#' # generating model
#' n=50
#' set.seed(8) # a seed of 2 also works
#' x1 = rnorm(n,0.22,2)
#' x7 = 0.5*x1 + rnorm(n,0,sd=2)
#' x6 = -0.75*x1 + rnorm(n,0,3)
#' x3 = -0.5-0.5*x6 + rnorm(n,0,2)
#' x9 = rnorm(n,0.6,3.5)
#' x4 = 0.5*x9 + rnorm(n,0,sd=3)
#' x2 = -0.5 + 0.5*x9 + rnorm(n,0,sd=2)
#' x5 = -0.5*x2+0.5*x3+0.5*x6-0.5*x9+rnorm(n,0,1.5)
#' x8 = x1 + x2 -2*x3 - 0.3*x4 + x5 - 1.6*x6 - 1*x7 + x9 +rnorm(n,0,0.5)
#' y = 0.6*x8 + rnorm(n,0,2)
#' artificialeg = round(data.frame(x1,x2,x3,x4,x5,x6,x7,x8,x9,y),1)
NULL




#' Forced Expiratory Volume
#'
#' This data set consists of 654 observations on youths aged 3 to 19 from 
#' East Boston recorded duing the middle to late 1970's. 
#' Forced expiratory volume (FEV), a measure of lung capacity, is the 
#' variable of interest. Age and height are two continuous predictors. 
#' Sex and smoke are two categorical predictors.
#'
#' @name fev
#' @format A data frame with 654 observations on 5 variables.
#' \describe{
#' \item{age}{Age (years)}
#' \item{fev}{Forced expiratory volume (liters).  Roughly the amount 
#'            of air an individual can exhale in the first second of 
#'            a forceful breath.}
#' \item{height}{Height (inches).}
#' \item{sex}{Female is 0. Male is 1.}
#' \item{smoke}{A binary variable indicating whether or not the 
#'              youth smokes. Nonsmoker is 0. Smoker is 1.}
#' }
#' @details Copies of this data set can also be found in the 
#'  \code{coneproj} and \code{tmle} packages.
#' @references 
#'  Tager, I. B., Weiss, S. T., Rosner, B., and Speizer, F. E. (1979). 
#'  Effect of parental cigarette smoking on pulmonary function in children. 
#'  \emph{American Journal of Epidemiology}, \bold{110}, 15-26.
#'  
#'  Rosner, B. (1999).
#'  \emph{Fundamentals of Biostatistics}, 5th Ed., Pacific Grove, CA: Duxbury.
#'   
#'   Kahn, M.J. (2005). An Exhalent Problem for Teaching Statistics.
#'   \emph{Journal of Statistics Education},  \bold{13}(2). 
#'    http://www.amstat.org/publications/jse/v13n2/datasets.kahn.html
#' @docType data
#' @keywords datasets
#' @usage data(fev)
#' @examples
#' data(fev)
#' full.mod = lm(fev~.,data=fev)
#' step(full.mod)
NULL
