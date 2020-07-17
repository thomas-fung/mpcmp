context("Test Update function")
library(mpcmp)

test_that("Test updating the mean regression formula", 
          { 
            data("takeoverbids")
            M.bids.full <- glm.cmp(numbids ~ leglrest + rearest + 
                                     finrest + whtknght + 
                                     bidprem + insthold + size + 
                                     sizesq + regulatn, data=takeoverbids)
            M.bids.null <- update(M.bids.full, .~.-whtknght)
            expect_equal(round(M.bids.null$residuals_deviance,3), 131.838)
            expect_equal(length(M.bids.full$coef), 10)
            expect_equal(length(M.bids.null$coef), 9)
          })


test_that("Test updating the dispersion regression formula",
          {
            data(sitophilus)
            M.sit.full <- glm.cmp(formula = ninsect ~ extract, 
                                  formula_nu = ~extract, data = sitophilus)
            M.sit.null1 <- update(M.sit.full, formula_nu. =  ~.-extract)
            expect_equal(round(M.sit.null1$residuals_deviance,3), 41.408)
            expect_equal(length(M.sit.full$coefficients), 8)
            expect_equal(length(M.sit.null1$coefficients), 5)
            expect_equal(length(M.sit.null1$nu), 40)
            M.sit.null2 <- update(M.sit.full, formula_nu. = NULL)
            expect_equal(length(M.sit.null2$nu), 1)
            data(attendance)
            M.attendance.fix.nu <- glm.cmp(daysabs~ gender+math+prog, 
                                           data=attendance)
            M.attendance.vary.nu <- update(M.attendance.fix.nu, 
                                           formula_nu. = ~ 1+math)
            expect_equal(length(M.attendance.fix.nu$nu),1)
            expect_equal(length(M.attendance.vary.nu$nu), 314)
          })


