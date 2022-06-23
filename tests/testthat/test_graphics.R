context("Test Graphics")
library(mpcmp)
data("attendance")
M.attendance <- glm.cmp(daysabs ~ gender + math + prog,
  data = attendance
)
test_that("Test gg_plot", {
  expect_length(gg_plot(M.attendance, which = 1:8), 8)
  expect_warning(gg_plot(M.attendance, which = 9),
                 "The acceptable ragne for option 'which' is 1:8.")
}
)
test_that("Test base_plot", {
  expect_warning(plot(M.attendance, which = 9),
                 "The acceptable ragne for option 'which' is 1:8.")
}
)
