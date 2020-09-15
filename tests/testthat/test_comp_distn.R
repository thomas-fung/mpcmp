context("Test CMP distribution related functions")
library(mpcmp)

test_that("Test pcomp", 
          { 
            expect_equal(round(pcomp(5, 4, 0.5),3), 0.743)
            expect_equal(pcomp(-1, 4, 0.5), 0)
            expect_equal(round(pcomp(5, 4, 1.5),3), 0.822)
            expect_error(pcomp(1.5), 'argument "mu" is missing, with no default')
            expect_error(pcomp(1.5, mu = 1.5, lambda = 1.5, nu = 0.9), 
                         "specify 'mu' or 'lambda' but not both")
          })


test_that("Test dcomp",
          {
            expect_equal(round(dcomp(5, 4, 0.5),3), 0.121)
            expect_warning(dcomp(5.5, 4, 0.5))
            expect_equal(dcomp(-1, 4, 0.5), 0)
            expect_error(dcomp(1.5), 'argument "mu" is missing, with no default')
            expect_error(dcomp(1.5, mu = 1.5, lambda = 1.5, nu = 0.9), 
                         "specify 'mu' or 'lambda' but not both")
          })

test_that("Test qcomp",
          {
            expect_equal(qcomp(0.2, 4, 0.5), 2)
            expect_equal(qcomp(0, 4, 0.5), 0)
            expect_warning(qcomp(1.1, 4, 0.5))
            expect_warning(qcomp(-1, 4, 0.5))
            expect_error(qcomp(1.5), 'argument "mu" is missing, with no default')
            expect_error(qcomp(1.5, mu = 1.5, lambda = 1.5, nu = 0.9), 
                         "specify 'mu' or 'lambda' but not both")
          })

test_that("Test rcomp",
          {
            expect_length(rcomp(20, 4, 0.5), 20)
            expect_error(rcomp(20), 'argument "mu" is missing, with no default')
            expect_error(rcomp(20, mu = 1.5, lambda = 1.5, nu = 0.9), 
                         "specify 'mu' or 'lambda' but not both")
          })
