context("Test Graphics")
library(mpcmp)
data("attendance")
M.attendance <- glm.cmp(daysabs~ gender+math+prog, 
                        data=attendance)
test_that("Test gg_plot", 
          { 
            expect_length(gg_plot(M.attendance, which = 1:8), 8)
            expect_output(gg_plot(M.attendance, which=9), 
                          NULL)
          })

