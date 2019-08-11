context("Testing mean approximation")

library(mpcmp)

test_that("Testing theoretical against numerical mean", { 
  t.mean1 <- sum(0:100*dcomp(0:100,2,0.2))
  t.mean2 <- sum(0:100*dcomp(0:100,40,0.2))
  t.mean3 <- sum(0:100*dcomp(0:100,2,4.5))
  t.mean4 <- sum(0:100*dcomp(0:100,40,4.5))
  t.mean5 <- sum(0:500*dcomp(0:500,150,0.8))
  expect_equal(round(abs(t.mean1-2)/2,3), 0)
  expect_equal(round(abs(t.mean2-40)/40,3), 0)
  expect_equal(round(abs(t.mean3-2)/2,3), 0)
  expect_equal(round(abs(t.mean4-40)/40,3), 0)
  expect_equal(round(abs(t.mean5-150)/150,3), 0)
})


