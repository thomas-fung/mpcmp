context("Test Graphics")
library(mpcmp)

test_that("Test gg_plot", 
          { data("attendance")
            M.attendance <- glm.cmp(daysabs~ gender+math+prog, 
                                    data=attendance)
            expect_length(gg_plot(M.attendance, which = 1:7), 7)
            expect_warning(gg_plot(M.attendance, which=8))
          })

