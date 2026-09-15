library(mpcmp)

# x3 is an exact linear combination of x1 and x2, so the mean (and,
# separately, the dispersion) design matrix used below is rank-deficient by
# construction -- this mirrors how glm() handles collinear designs, and is
# what fit_glm_cmp_const_nu()/fit_glm_cmp_vary_nu() are meant to replicate.
set.seed(4821)
n <- 100
x1 <- rnorm(n)
x2 <- rnorm(n)
x3 <- x1 + 2 * x2
mu <- exp(0.5 + 0.3 * x1 - 0.2 * x2)
y <- rpois(n, mu)
dat_rd <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3)

M.rd <- suppressWarnings(glm.cmp(y ~ x1 + x2 + x3, data = dat_rd))
M.rd.null <- glm.cmp(y ~ x1 + x2, data = dat_rd)
M.rd.vary <- suppressWarnings(
  glm.cmp(y ~ x1 + x2, formula_nu = ~ x1 + x2 + x3, data = dat_rd)
)

test_that("glm.cmp warns and reports NA for an aliased mean coefficient", {
  expect_warning(
    glm.cmp(y ~ x1 + x2 + x3, data = dat_rd),
    "rank-deficient"
  )
  expect_equal(M.rd$rank, 3)
  expect_true(is.na(coef(M.rd)["x3"]))
  expect_false(anyNA(coef(M.rd)[c("(Intercept)", "x1", "x2")]))
})

test_that("glm.cmp warns and reports NA for an aliased dispersion coefficient", {
  expect_warning(
    glm.cmp(y ~ x1 + x2, formula_nu = ~ x1 + x2 + x3, data = dat_rd),
    "rank-deficient"
  )
  expect_equal(M.rd.vary$rank_nu, 3)
  expect_true(is.na(M.rd.vary$coefficients_gamma["x3"]))
  expect_false(anyNA(M.rd.vary$coefficients_gamma[c("(Intercept)", "x1", "x2")]))
})

test_that("model.matrix.cmp returns the full design matrix, matching length(coef)", {
  mm <- model.matrix(M.rd)
  expect_equal(ncol(mm), length(coef(M.rd)))
  expect_equal(colnames(mm), names(coef(M.rd)))
  expect_equal(nrow(mm), nobs(M.rd))

  # a full-rank fit is unaffected
  mm_full <- model.matrix(M.rd.null)
  expect_equal(ncol(mm_full), length(coef(M.rd.null)))
  expect_false(anyNA(coef(M.rd.null)))

  mm_vary <- model.matrix(M.rd.vary)
  expect_equal(ncol(mm_vary$x), length(M.rd.vary$coefficients_beta))
  expect_equal(ncol(mm_vary$s), length(M.rd.vary$coefficients_gamma))
})

test_that("cmplrtest uses numeric rank rather than raw coefficient length", {
  # Because x3 adds nothing over x1 + x2, the rank-deficient fit is
  # statistically equivalent to the nested (full-rank) null model: the LR
  # test should report 0 additional estimated parameters and (up to
  # numerical error) a zero test statistic. Before the rank fix,
  # length(coefficients) would have (incorrectly) reported df = 1, since
  # the aliased NA coefficient was counted as an extra parameter.
  lrt <- cmplrtest(M.rd, M.rd.null)
  expect_s3_class(lrt, "htest")
  expect_equal(unname(lrt$parameter), 0)
  expect_equal(unname(lrt$statistic), 0, tolerance = 1e-6)
})

test_that("confint.cmp reports NA only for the aliased coefficient", {
  ci <- confint(M.rd)
  expect_true(is.matrix(ci))
  expect_true(all(is.na(ci["x3", ])))
  expect_false(anyNA(ci[c("(Intercept)", "x1", "x2"), ]))

  ci_full <- confint(M.rd.null)
  expect_false(anyNA(ci_full))
})
