library(mpcmp)

data(sitophilus)

# Regression test: fit_glm_cmp_vary_nu() used to error with
# "object 'lambdaold' not found" when 'gammastart' was supplied, because the
# initial comp_lambdas() call referenced a variable that only exists inside
# the iteration loop.
test_that("glm.cmp with formula_nu runs when gammastart is supplied", {
  expect_no_error(
    M <- glm.cmp(
      ninsect ~ extract,
      formula_nu = ~extract,
      gammastart = rep(0, 4),
      data = sitophilus
    )
  )
})

M.sit <- glm.cmp(
  ninsect ~ extract,
  formula_nu = ~extract,
  gammastart = rep(0, 4),
  data = sitophilus
)
M.sit.nostart <- glm.cmp(
  ninsect ~ extract,
  formula_nu = ~extract,
  data = sitophilus
)

test_that("gammastart fit returns a valid cmp object with full-length coefs", {
  expect_s3_class(M.sit, "cmp")
  expect_length(M.sit$coefficients_beta, 4)
  expect_length(M.sit$coefficients_gamma, 4)
  expect_false(anyNA(M.sit$coefficients_gamma))
})

test_that("gammastart fit converges to the same solution as the default start", {
  expect_equal(M.sit$maxl, M.sit.nostart$maxl, tolerance = 1e-4)
  expect_equal(coef(M.sit), coef(M.sit.nostart), tolerance = 1e-3)
})

test_that("gammastart of the wrong length errors informatively", {
  expect_error(
    glm.cmp(
      ninsect ~ extract,
      formula_nu = ~extract,
      gammastart = rep(0, 2),
      data = sitophilus
    ),
    "gammastart"
  )
})
