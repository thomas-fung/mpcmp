context("Test estimated coefficients")
library(mpcmp)

test_that("Test the estimated coefficients from the attendance dataset", 
          { 
            data("attendance")
            M.attendance <- glm.cmp(daysabs~ gender+math+prog, data=attendance)
            expect_equal(unname(round(M.attendance$nu,4)), 0.0202)
            expect_equal(unname(round(coef(M.attendance)[1],3)), 2.715)
          })
