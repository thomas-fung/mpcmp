context("Test Summarise and Extract info from CMP object")
library(mpcmp)
data("attendance")
M.attendance <- glm.cmp(daysabs~ gender+math+prog, 
                        data=attendance)

test_that("Test residuals", 
          { 
            expect_equal(-0.250, 
                         unname(round(residuals(M.attendance)[1],3)))
            expect_equal(unname(round(residuals(M.attendance, 
                                         type = "pearson")[1],3)), -0.232)
            expect_equal(unname(round(residuals(M.attendance, 
                                         type = "response")[1],3)), -1.249)
          })

test_that("Test nobs",
          {
            expect_equal(nobs(M.attendance), 314)
          })

test_that("Test fitted",
          {
            expect_length(fitted(M.attendance), 314)
          })

test_that("Test Model frame",
          {
            expect_equal(NROW(model.frame(M.attendance)), 314)
            data("sitophilus")
            data(sitophilus)
            M.sit <- glm.cmp(formula = ninsect ~ extract, 
                             formula_nu = ~extract, data = sitophilus)
            expect_length(model.frame(M.sit),2)
          })

test_that("Test coefficients",
          {
            expect_length(coefficients(M.attendance),5)
          })

test_that("Test predict",
          {
            expect_length(predict(M.attendance), 314)
            expect_length(predict(M.attendance, se.fit = TRUE),2)
            expect_length(predict(M.attendance, type = "response"), 314)
            expect_length(predict(M.attendance, type = "response", 
                                  se.fit = TRUE),2)
            expect_length(predict(M.attendance, 
                                  newdata = attendance[2,]),1)
            expect_length(predict(M.attendance, 
                                  newdata = attendance[1:2,],
                                  se.fit = TRUE), 2)
          })
