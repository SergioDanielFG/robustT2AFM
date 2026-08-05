make_cal <- function(w = c(B1 = 0.10, B2 = 0.40, B3 = 0.25, B4 = 0.25)) {
  list(weights = w, mu_r = c(0, 0), Sw = diag(2))
}

test_that("plot_afm_weights returns a ggplot sorted from lowest weight", {
  p <- plot_afm_weights(make_cal())

  expect_s3_class(p, "ggplot")

  # One bar per batch, ascending by weight.
  expect_equal(nrow(p$data), 4L)
  expect_equal(p$data$Weight, c(0.10, 0.25, 0.25, 0.40))
  expect_equal(as.character(p$data$Batch)[1], "B1")

  # Factor levels are reversed so the smallest weight is drawn at the top:
  # the first level is the one ggplot places at the bottom.
  expect_equal(levels(p$data$Batch)[1], "B2")
  expect_equal(levels(p$data$Batch)[4], "B1")
})

test_that("plot_afm_weights marks the uniform weight 1/K", {
  p <- plot_afm_weights(make_cal())

  # The dashed reference line sits exactly at 1/K.
  vline <- Filter(function(l) inherits(l$geom, "GeomVline"), p$layers)
  expect_length(vline, 1L)
  expect_equal(vline[[1]]$data$xintercept, 1 / 4)

  # And it is labelled with the value.
  expect_true(any(grepl("1/K = 0.2500", vapply(p$layers, function(l) {
    lbl <- l$aes_params$label
    if (is.null(lbl)) "" else as.character(lbl)
  }, character(1)), fixed = TRUE)))
})

test_that("plot_afm_weights highlights only what the caller asks for", {
  # Default: every bar neutral, no batch singled out by the package.
  p0 <- plot_afm_weights(make_cal())
  expect_false(any(p0$data$Highlight))
  expect_equal(length(unique(p0$data$Fill)), 1L)

  # highlight_lowest = 2 marks the two smallest weights and nothing else.
  p2 <- plot_afm_weights(make_cal(), highlight_lowest = 2)
  expect_equal(sum(p2$data$Highlight), 2L)
  expect_true(all(p2$data$Highlight[1:2]))
  expect_equal(as.character(p2$data$Batch)[p2$data$Highlight],
               c("B1", "B3"))
})

test_that("plot_afm_weights works on a real calibration object", {
  sim  <- simulate_batch_process(K1 = 12, K2 = 0, I = 20, J = 4,
                                 seed = 20260417)
  cal  <- calibrate_afm_mcd(sim, paste0("Var", 1:4))
  p    <- plot_afm_weights(cal)

  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), length(cal$weights))
  expect_false(is.unsorted(p$data$Weight))
})

test_that("plot_afm_weights writes PNG and PDF only when save_path is given", {
  stem <- file.path(tempdir(), "test_afm_weights")
  on.exit(unlink(paste0(stem, c(".png", ".pdf"))), add = TRUE)

  # Default writes nothing to disk (CRAN policy).
  invisible(plot_afm_weights(make_cal()))
  expect_false(file.exists(paste0(stem, ".png")))

  invisible(plot_afm_weights(make_cal(), save_path = stem))
  expect_true(file.exists(paste0(stem, ".png")))
  expect_true(file.exists(paste0(stem, ".pdf")))
})

test_that("plot_afm_weights validates its inputs", {
  expect_error(plot_afm_weights("nope"), "must be a list from calibrate_afm_mcd")
  expect_error(plot_afm_weights(list(mu_r = 1)),
               "must be a list from calibrate_afm_mcd")
  expect_error(plot_afm_weights(list(weights = c(B1 = 1))),
               "at least 2 finite")
  expect_error(plot_afm_weights(list(weights = c(0.5, 0.5))),
               "must be named with the batch identifiers")
  expect_error(plot_afm_weights(make_cal(), title = c("a", "b")),
               "'title' must be")
  expect_error(plot_afm_weights(make_cal(), highlight_lowest = -1),
               "non-negative whole number")
  expect_error(plot_afm_weights(make_cal(), highlight_lowest = 1.5),
               "non-negative whole number")
  expect_error(plot_afm_weights(make_cal(), highlight_lowest = 99),
               "exceeds the number of")
  expect_error(plot_afm_weights(make_cal(), save_path = ""), "save_path")
})
