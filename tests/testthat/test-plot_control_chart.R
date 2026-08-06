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

test_that("plot_control_chart is quiet by default and verbose on request", {
  mon <- data.frame(
    Batch = paste0("B", 1:5),
    I     = rep(20L, 5),
    T2    = c(3, 5, 25, 8, 30)
  )
  n_text <- function(p) {
    sum(vapply(p$layers, function(l) inherits(l$geom, "GeomText"), logical(1)))
  }

  # Default: no per-point value labels. Three text layers remain: the OOC
  # count, the UCL label and the identifiers of the OOC batches.
  p0 <- plot_control_chart(mon, UCL = 20)
  expect_equal(n_text(p0), 3L)

  # show_values = TRUE adds a fourth layer, one label per batch.
  p1 <- plot_control_chart(mon, UCL = 20, show_values = TRUE)
  expect_equal(n_text(p1), 4L)

  # The legend sits below the panel, as in plot_method_comparison.
  expect_equal(p0$theme$legend.position, "bottom")

  # Dropping the value labels frees head-room: the default y range is lower.
  expect_lt(max(p0$scales$get_scales("y")$limits),
            max(p1$scales$get_scales("y")$limits))
})

test_that("out-of-control labels are staggered only when batches are close", {
  # B1..B3 consecutive and out of control, then a gap, then B8.
  mon <- data.frame(
    Batch = paste0("B", 1:8),
    I     = rep(20L, 8),
    T2    = c(25, 26, 27, 2, 3, 2, 4, 30)
  )
  p <- plot_control_chart(mon, UCL = 20)

  lab <- p$layers[[which(vapply(p$layers, function(l) {
    "vjust" %in% names(l$data)
  }, logical(1)))[1]]]$data
  lab <- lab[order(lab$BatchIdx), ]

  expect_equal(as.character(lab$Batch), c("B1", "B2", "B3", "B8"))
  # Alternating inside the run, back to the base level after the gap.
  expect_equal(lab$vjust[1], lab$vjust[3])
  expect_false(lab$vjust[1] == lab$vjust[2])
  expect_equal(lab$vjust[4], lab$vjust[1])
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
