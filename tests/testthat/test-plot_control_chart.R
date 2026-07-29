test_that("plot_control_chart returns a ggplot object", {
  mon <- data.frame(
    Batch = paste0("B", 1:5),
    I     = rep(20L, 5),
    T2    = c(3, 5, 25, 8, 30)   # two batches above UCL = 20
  )
  p <- plot_control_chart(mon, UCL = 20,
                          method_label = "Test chart",
                          alpha = 0.001)
  expect_s3_class(p, "ggplot")
})

test_that("plot_control_chart writes PNG and PDF when save_path is given", {
  mon <- data.frame(
    Batch = paste0("B", 1:5),
    I     = rep(20L, 5),
    T2    = c(3, 5, 25, 8, 30)
  )
  stem <- file.path(tempdir(), "test_control_chart")
  on.exit(unlink(paste0(stem, c(".png", ".pdf"))), add = TRUE)

  invisible(plot_control_chart(mon, UCL = 20,
                               method_label = "Test chart",
                               alpha = 0.001,
                               save_path = stem))
  expect_true(file.exists(paste0(stem, ".png")))
  expect_true(file.exists(paste0(stem, ".pdf")))
})

test_that("plot_control_chart validates its inputs", {
  mon <- data.frame(Batch = "B1", T2 = 1)
  expect_error(plot_control_chart("nope", UCL = 20), "must be a data frame")
  expect_error(plot_control_chart(mon,   UCL = -1),  "positive")
  expect_error(plot_control_chart(mon,   UCL = 20, alpha = 2),
               "alpha' must be")
  expect_error(plot_control_chart(mon,   UCL = 20, save_path = ""),
               "save_path")
})
