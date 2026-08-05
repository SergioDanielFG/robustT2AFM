# Synthetic pair of monitoring results: batch B3 is signalled by the robust
# chart only, B5 by the classical chart only, B1 sits inside the near-miss
# band of the classical limit (0.85 * 20 = 17).
make_pieces <- function() {
  list(
    mon_r = data.frame(Batch = paste0("B", 1:5), I = rep(20L, 5),
                       T2 = c(10, 5, 30, 2, 8),
                       stringsAsFactors = FALSE),
    mon_c = data.frame(Batch = paste0("B", 1:5), I = rep(20L, 5),
                       T2 = c(17, 4, 12, 3, 25),
                       stringsAsFactors = FALSE),
    ucl_r = 20, ucl_c = 20
  )
}

# The layer carrying one row per batch per method.
points_layer <- function(p) {
  i <- which(vapply(p$layers, function(l) "flagged" %in% names(l$data),
                    logical(1)))[1]
  p$layers[[i]]$data
}

make_study <- function(compare = TRUE) {
  sim <- simulate_batch_process(
    K1 = 12, K2 = 8, I = 20, J = 4, rho = 0.6,
    outlier_batches_F1 = 3, outlier_rate = 0.20, outlier_shift = 4,
    prop_ooc_F2 = 0.5, shift_ooc = 1.0, seed = 20260425
  )
  run_afm_mcd(subset(sim, Phase == "Phase 1"), subset(sim, Phase == "Phase 2"),
              paste0("Var", 1:4), plot = FALSE, compare_classical = compare)
}

test_that("plot_method_comparison accepts both input forms", {
  pz <- make_pieces()
  p1 <- plot_method_comparison(pz$mon_r, pz$mon_c, pz$ucl_r, pz$ucl_c)
  expect_s3_class(p1, "ggplot")

  study <- make_study()
  p2 <- plot_method_comparison(study)
  expect_s3_class(p2, "ggplot")
  # One point per batch per method: the two panels cover the same batches.
  pts <- points_layer(p2)
  expect_equal(nrow(pts), 2L * nrow(study$monitoring))
  expect_equal(levels(pts$Method),
               c("AFM-MCD (robust)", "Classical Hotelling"))

  # A study without the classical baseline says what to re-run.
  expect_error(plot_method_comparison(make_study(compare = FALSE)),
               "compare_classical = TRUE")
})

test_that("the default scale is the ratio to each method's own limit", {
  pz <- make_pieces()
  p  <- plot_method_comparison(pz$mon_r, pz$mon_c, pz$ucl_r, pz$ucl_c)

  pts <- points_layer(p)
  expect_equal(pts$y, pts$T2 / pts$UCL)

  # Raw statistics remain available.
  p_t2 <- plot_method_comparison(pz$mon_r, pz$mon_c, pz$ucl_r, pz$ucl_c,
                                 scale = "T2")
  expect_equal(points_layer(p_t2)$y, points_layer(p_t2)$T2)

  # Only the raw-scale caption warns that heights are not comparable.
  expect_true(grepl("NOT\ncomparable|NOT comparable", p_t2$labels$caption))
  expect_false(grepl("NOT comparable", p$labels$caption))
  expect_true(grepl("share one scale", p$labels$caption))
})

test_that("disagreeing batches are ringed in the panel that stays silent", {
  pz <- make_pieces()
  p  <- plot_method_comparison(pz$mon_r, pz$mon_c, pz$ucl_r, pz$ucl_c)

  ring <- p$layers[[which(vapply(p$layers, function(l) {
    isTRUE(l$aes_params$shape == 21)
  }, logical(1)))[1]]]$data

  # B3: robust only -> ringed in the classical panel.
  # B5: classical only -> ringed in the robust panel.
  expect_setequal(ring$Batch, c("B3", "B5"))
  expect_equal(ring$Method[ring$Batch == "B3"],
               factor("Classical Hotelling",
                      levels = c("AFM-MCD (robust)", "Classical Hotelling")))
  expect_equal(ring$Method[ring$Batch == "B5"],
               factor("AFM-MCD (robust)",
                      levels = c("AFM-MCD (robust)", "Classical Hotelling")))
})

test_that("without faulty the figure talks about disagreement, never misses", {
  pz  <- make_pieces()
  p   <- plot_method_comparison(pz$mon_r, pz$mon_c, pz$ucl_r, pz$ucl_c)
  txt <- paste(p$labels$caption,
               paste(unlist(lapply(p$layers, function(l) l$data$label)),
                     collapse = " "))

  expect_false(grepl("missed", txt, ignore.case = TRUE))
  expect_false(grepl("faulty", txt, ignore.case = TRUE))
  expect_true(grepl("signalled by one method only", txt))
  # Robust signals B3 only; classical signals B5 only.
  expect_true(grepl("signalled 1 of 5 batches", txt))
})

test_that("with faulty the figure counts detections and may say missed", {
  pz  <- make_pieces()
  p   <- plot_method_comparison(pz$mon_r, pz$mon_c, pz$ucl_r, pz$ucl_c,
                                faulty = c("B3", "B5"))
  txt <- paste(p$labels$caption,
               paste(unlist(lapply(p$layers, function(l) l$data$label)),
                     collapse = " "))

  expect_true(grepl("detected 1 of 2 faulty", txt))
  expect_true(grepl("false alarm", txt))
  expect_true(grepl("missed by the classical chart", txt))

  # A typo in 'faulty' is caught rather than silently ignored.
  expect_error(plot_method_comparison(pz$mon_r, pz$mon_c, pz$ucl_r, pz$ucl_c,
                                      faulty = "B99"),
               "not among the monitored batches")
})

test_that("the near-miss band counts batches and never claims a near detection", {
  pz <- make_pieces()
  p  <- plot_method_comparison(pz$mon_r, pz$mon_c, pz$ucl_r, pz$ucl_c)

  # B1 has classical T2 = 17, inside [0.8*20, 20].
  expect_true(grepl("within 20% below the classical limit \\(1 batch",
                    p$labels$caption))
  expect_true(grepl("does not mean\nthey were almost detected|does not mean they were almost detected",
                    p$labels$caption))

  # The band is a single rect confined to the classical panel.
  band <- p$layers[[1]]$data
  expect_equal(as.character(band$Method), "Classical Hotelling")
  expect_equal(band$ymin, 0.8)
  expect_equal(band$ymax, 1)

  # near_miss = NULL removes both the band and its caption line.
  p0 <- plot_method_comparison(pz$mon_r, pz$mon_c, pz$ucl_r, pz$ucl_c,
                               near_miss = NULL)
  expect_false(grepl("Shaded band", p0$labels$caption))
})

test_that("plot_method_comparison writes files only when asked", {
  pz   <- make_pieces()
  stem <- file.path(tempdir(), "test_method_comparison")
  on.exit(unlink(paste0(stem, c(".png", ".pdf"))), add = TRUE)

  invisible(plot_method_comparison(pz$mon_r, pz$mon_c, pz$ucl_r, pz$ucl_c))
  expect_false(file.exists(paste0(stem, ".png")))

  invisible(plot_method_comparison(pz$mon_r, pz$mon_c, pz$ucl_r, pz$ucl_c,
                                   save_path = stem))
  expect_true(file.exists(paste0(stem, ".png")))
  expect_true(file.exists(paste0(stem, ".pdf")))
})

test_that("plot_method_comparison validates its inputs", {
  pz <- make_pieces()

  expect_error(plot_method_comparison("nope"), "must be an afm_mcd_study")
  expect_error(plot_method_comparison(pz$mon_r), "'classical' must be")
  expect_error(plot_method_comparison(pz$mon_r, pz$mon_c, -1, 20),
               "'UCL_robust' must be")
  expect_error(plot_method_comparison(pz$mon_r, pz$mon_c, 20, 0),
               "'UCL_classical' must be")
  expect_error(plot_method_comparison(pz$mon_r, pz$mon_c[1:3, ], 20, 20),
               "different batches")
  expect_error(plot_method_comparison(pz$mon_r, pz$mon_c, 20, 20,
                                      near_miss = 1.5),
               "'near_miss' must be")
  expect_error(plot_method_comparison(pz$mon_r, pz$mon_c, 20, 20,
                                      labels = c("a", "a")),
               "two distinct character strings")
  expect_error(plot_method_comparison(pz$mon_r, pz$mon_c, 20, 20,
                                      save_path = ""),
               "save_path")
  expect_error(plot_method_comparison(pz$mon_r, pz$mon_c, 20, 20,
                                      scale = "percent"),
               "'arg' should be one of")
})
