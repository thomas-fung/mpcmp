context("Test Inference")
library(mpcmp)
data("takeoverbids")
data("sitophilus")
M.sit.full <- glm.cmp(formula = ninsect ~ extract, 
                      formula_nu = ~extract, data = sitophilus)
M.bids.full <- glm.cmp(numbids ~ leglrest + rearest + 
                         finrest + whtknght + 
                         bidprem + insthold + size + 
                         sizesq + regulatn, data=takeoverbids)
M.attendance <- glm.cmp(daysabs~ gender+math+prog, data=attendance)
M.bids.null <- update(M.bids.full, .~.-whtknght)

test_that("Test updating the mean regression formula", 
          { 
            expect_equal(round(M.bids.null$residual_deviance,3), 131.838)
            expect_equal(length(M.bids.full$coef), 10)
            expect_equal(length(M.bids.null$coef), 9)
          })

test_that("Test cmplrtest function",
          {
            expect_equal(capture_output_lines(cmplrtest(M.bids.full, M.bids.null))[5],
                         "P-value:  0.000214 ")
          }
          )

test_that("Test LRTnu", 
          {
            expect_equal(capture_output_lines(LRTnu(M.attendance))[8], 
                         "P-value: < 2e-16")
          })

test_that("Test updating the dispersion regression formula",
          {
            M.sit.null1 <- update(M.sit.full, formula_nu. =  ~.-extract)
            expect_equal(round(M.sit.null1$residual_deviance,3), 41.408)
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

test_that("Test the confint function",
          {
            expect_is(confint.cmp(M.attendance), class = "matrix")
            expect_equal(colnames(confint.cmp(M.attendance)), 
                         c("2.5 %", "97.5 %"))
            expect_equal(colnames(confint.cmp(M.attendance, parm = "math", level = 0.9)), 
                         c("5 %", "95 %"))
            expect_length(confint.cmp(M.attendance, parm = "math", level = 0.9), 2)
          })
