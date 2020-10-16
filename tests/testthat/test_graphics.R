context("Test Graphics")
library(mpcmp)
data("attendance")
M.attendance <- glm.cmp(daysabs~ gender+math+prog, 
                        data=attendance)
test_that("Test gg_plot", 
          { 
            expect_length(gg_plot(M.attendance, which = 1:8), 8)
            expect_output(gg_plot(M.attendance, which=9), 
                          NULL)
            disp_ggplot1 <- gg_plot(M.attendance, which = 1)
            disp_ggplot2 <- gg_plot(M.attendance, which = 2)
            disp_ggplot3 <- gg_plot(M.attendance, which = 3)
            set.seed(200253046)
            disp_ggplot4 <- gg_plot(M.attendance, which = 4)
            set.seed(200253046)
            disp_ggplot5 <- gg_plot(M.attendance, which = 5)
            disp_ggplot6 <- gg_plot(M.attendance, which = 6)
            disp_ggplot7 <- gg_plot(M.attendance, which = 7)
            disp_ggplot8 <- gg_plot(M.attendance, which = 8)
            vdiffr::expect_doppelganger("ggplot2 diagnostic 1",
                                        disp_ggplot1)
            vdiffr::expect_doppelganger("ggplot2 diagnostic 2",
                                        disp_ggplot2)
            vdiffr::expect_doppelganger("ggplot2 diagnostic 3",
                                        disp_ggplot3)
            vdiffr::expect_doppelganger("ggplot2 diagnostic 4",
                                        disp_ggplot4)
            vdiffr::expect_doppelganger("ggplot2 diagnostic 5",
                                        disp_ggplot5)
            vdiffr::expect_doppelganger("ggplot2 diagnostic 6",
                                        disp_ggplot6)
            vdiffr::expect_doppelganger("ggplot2 diagnostic 7",
                                        disp_ggplot7)
            vdiffr::expect_doppelganger("ggplot2 diagnostic 8",
                                        disp_ggplot8)
          })

test_that("Test plot", {
  disp_plot1 <- function() plot.cmp(M.attendance, which = 1)
  disp_plot2 <- function() plot.cmp(M.attendance, which = 2)
  disp_plot3 <- function() plot.cmp(M.attendance, which = 3)
  disp_plot4 <- function() {
    set.seed(200253046)
    plot.cmp(M.attendance, which = 4)
    }
  disp_plot5 <- function() {
    set.seed(200253046)
    plot.cmp(M.attendance, which = 5)
  }
  disp_plot6 <- function() plot.cmp(M.attendance, which = 6)
  disp_plot7 <- function() plot.cmp(M.attendance, which = 7)
  disp_plot8 <- function() plot.cmp(M.attendance, which = 8)
  vdiffr::expect_doppelganger("plot diagnostic1", disp_plot1)
  vdiffr::expect_doppelganger("plot diagnostic2", disp_plot2)
  vdiffr::expect_doppelganger("plot diagnostic3", disp_plot3)
  vdiffr::expect_doppelganger("plot diagnostic4", disp_plot4)
  vdiffr::expect_doppelganger("plot diagnostic5", disp_plot5)
  vdiffr::expect_doppelganger("plot diagnostic6", disp_plot6)
  vdiffr::expect_doppelganger("plot diagnostic7", disp_plot7)
  vdiffr::expect_doppelganger("plot diagnostic8", disp_plot8)
})