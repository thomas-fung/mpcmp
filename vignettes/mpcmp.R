## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library("mpcmp")
data("attendance", package="mpcmp")
M.attendance <- glm.cmp(daysabs~ gender+math+prog, data=attendance)
M.attendance
summary(M.attendance)

## ---- fig.show='hold', fig.cap='test'------------------------------------
plot(1:10)
plot(10:1)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(mtcars, 10))

