context("Test adding an offset term")
library(mpcmp)

test_that("Test adding an offset term", 
          {
            data("cottonbolls")
            set.seed(1)
            cottonbolls$n_samples <- rpois(nrow(cottonbolls), 5)
            M.bolls <- glm.cmp(nc~ 1+stages:def+stages:def2, offset = log(n_samples),data = cottonbolls)
            expect_equal(M.bolls$offset, log(cottonbolls$n_samples))
          })
