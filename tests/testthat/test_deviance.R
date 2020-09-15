context("Deviance output")
library(mpcmp)

test_that("Testing attendance example", 
          { 
            data(attendance)
            M.attendance <- glm.cmp(daysabs~ gender+math+prog, data=attendance)
            expect_equal(round(AIC.cmp(M.attendance),3),1739.026)
            expect_equal(round(M.attendance$null_deviance,3), 455.833)
            expect_equal(round(M.attendance$residual_deviance,3), 377.441)
          })


test_that("Testing takeover example",
          {
            data(takeoverbids)
            M.bids <- glm.cmp(numbids ~ leglrest + rearest + finrest + whtknght 
                              + bidprem + insthold + size + sizesq + regulatn, 
                              data=takeoverbids)
            expect_equal(round(AIC.cmp(M.bids),3),  382.175)
            expect_equal(round(M.bids$null_deviance,2), 182.39)
            expect_equal(round(M.bids$residual_deviance,1), 131.2)
          })

test_that("Testing cottonbolls example",
          { data(cottonbolls)
            M.bolls <- glm.cmp(nc~ 1+stages:def+stages:def2, data= cottonbolls)
            expect_equal(round(AIC.cmp(M.bolls),3), 440.823)
            expect_equal(round(M.bolls$null_deviance,2), 345.94)
            expect_equal(round(M.bolls$residual_deviance,1), 125.3)
          })

test_that("Testing fish example",
          { data(fish)
            M.fish <- glm.cmp(species~ 1+log(area), data=fish)
            expect_equal(round(AIC.cmp(M.fish),3),  638.853)
            expect_equal(round(M.fish$null_deviance,2), 101.67)
            expect_equal(round(M.fish$residual_deviance,1), 59.5)
          })

test_that("Testing sitophilus example",
          { data(sitophilus)
            M.sit <- glm.cmp(formula = ninsect ~ extract, 
                      formula_nu = ~extract, data = sitophilus)
            expect_equal(round(AIC.cmp(M.sit),3), 260.828)
            expect_equal(round(M.sit$null_deviance,2), 257.21)
            expect_equal(round(M.sit$residual_deviance,2), 41.68)
            })