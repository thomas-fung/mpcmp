context("Test Summarise and Extract info from CMP object")
library(mpcmp)
data("attendance")
M.attendance <- glm.cmp(daysabs~ gender+math+prog, 
                        data=attendance)
M.sit <- glm.cmp(formula = ninsect ~ extract, formula_nu = ~extract, data = sitophilus)

test_that("Test residuals", 
          { 
            expect_equal(-0.264, 
                         unname(round(residuals(M.attendance)[1],3)))
            expect_equal(unname(round(residuals(M.attendance, 
                                         type = "pearson")[1],3)), -0.244)
            expect_equal(unname(round(residuals(M.attendance, 
                                         type = "response")[1],3)), -1.345)
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

test_that("Test print",
          {
            expect_equal(capture_output_lines(print(M.attendance))[8], 
                         "Dispersion (nu): 0.0202")
          })


test_that("Test sumamry",
          {
            expect_equal(capture_output_lines(print(summary(M.attendance)))[24], 
                         "AIC: 1739.026 ")
            expect_is(summary(M.attendance)$coefficients, "matrix")
            expect_vector(as.vector(summary(M.attendance)$coefficients), 
                          ptype = numeric(), size = 20)
            expect_is(summary(M.sit)$coefficients, "matrix")
            expect_is(summary(M.sit)$coef.table_beta, "matrix")
            expect_is(summary(M.sit)$coef.table_gamma, "matrix")
            expect_equal(capture_output_lines(print(summary(M.sit)))[20],
                         "extractLeaf    -0.3831  0.6509  -0.589    0.556")
          })

test_that("Test rstandard",
         {
           expect_equal(round(unname(rstandard.cmp(M.attendance)[1]), 5), 
                        -0.26597)
           expect_equal(round(unname(rstandard.cmp(M.attendance, type = "pearson")[1]), 5), 
                        -0.2451)
         })

test_that("Test influence",
         {
           infl <- influence.cmp(M.attendance)
           expect_is(infl, class = "list")
           expect_equal(unname(infl$h[1]), 0.01152283)
           expect_equal(unname(round(infl$dev_res[1],4)), -0.2644)
           expect_equal(unname(round(infl$pear_res[1],4)), -0.2437)
         })

test_that("Test hatvalues",
          {
            expect_equal(unname(hatvalues.cmp(M.attendance)[1]),  0.01152283)
            expect_length(hatvalues.cmp(M.attendance), 314)
          })

test_that("Test cooks.distance",
          {
            expect_equal(unname(cooks.distance.cmp(M.attendance)[1]), 0.0001400528)
            expect_length(cooks.distance.cmp(M.attendance), 314)
          })

test_that("Test vcov",
          {
            expect_is(vcov(M.attendance), class = "matrix")
            expect_is(vcov(M.sit), class = "list")
          })

test_that("Test broom",
          {
            expect_named(tidy(M.sit), 
                         expected = c("parameter", "term", "estimate", 
                                      "std.error", "statistic", 
                                      "p.value"))
            expect_named(tidy(M.attendance), 
                         expected = c("term", "estimate", 
                                      "std.error", "statistic", 
                                      "p.value"))
            expect_length(glance(M.attendance), 8)
            expect_is(glance(M.attendance), class = "tbl_df")
            expect_length(glance(M.sit), 8)
            expect_is(glance(M.sit), class = "tbl_df")
            expect_is(augment(M.attendance), class = "tbl_df")
            expect_length(augment(M.attendance), 9)
            expect_length(augment(M.sit), 9)
          })

