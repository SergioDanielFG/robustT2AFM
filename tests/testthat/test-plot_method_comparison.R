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
n_geom <- function(p, cls) {
  sum(vapply(p$layers, function(l) inherits(l$geom, cls), logical(1)))
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
  expect_equal(nrow(points_layer(p2)), 2L * nrow(study$monitoring))
  expect_equal(levels(points_layer(p2)$Method),
               c("AFM-MCD (robust)", "Classical Hotelling"))

  # A study without the classical baseline says what to re-run.
  expect_error(plot_method_comparison(make_study(compare = FALSE)),
               "compare_classical = TRUE")
})

test_that("the default figure is bare: no caption and no marks", {
  pz <- make_pieces()
  p  <- plot_method_comparison(pz$mon_r, pz$mon_c, pz$ucl_r, pz$ucl_c)

  # No caption at all on the default scale.
  expect_null(p$labels$caption)

  # Exactly three layers: the limit line, the points, the limit label.
  expect_equal(length(p$layers), 3L)
  expect_equal(n_geom(p, "GeomText"), 1L)     # only the "own limit" label
  expect_equal(n_geom(p, "GeomVline"), 0L)    # no guides
  expect_equal(n_geom(p, "GeomRect"), 0L)     # no near-miss band

  # Two panels side by side.
  expect_equal(p$facet$params$ncol, 2L)
})

test_that("diagnostics adds guides, rings and counts, and nothing else does", {
  pz <- make_pieces()
  p  <- plot_method_comparison(pz$mon_r, pz$mon_c, pz$ucl_r, pz$ucl_c,
                               diagnostics = TRUE)

  expect_equal(n_geom(p, "GeomVline"), 1L)
  expect_equal(n_geom(p, "GeomText"), 2L)     # limit label + panel counts
  expect_null(p$labels$caption)               # still no caption

  ring <- p$layers[[which(vapply(p$layers, function(l) {
    isTRUE(l$aes_params$shape == 21)
  }, logical(1)))[1]]]$data

  # B3: robust only -> ringed in the classical panel.
  # B5: classical only -> ringed in the robust panel.
  expect_setequal(ring$Batch, c("B3", "B5"))
  expect_equal(as.character(ring$Method[ring$Batch == "B3"]),
               "Classical Hotelling")
  expect_equal(as.character(ring$Method[ring$Batch == "B5"]),
               "AFM-MCD (robust)")

  counts <- unlist(lapply(p$layers, function(l) l$data$label))
  expect_true(any(grepl("signalled 1 of 5", counts)))
  expect_false(any(grepl("faulty", counts)))  # no truth was supplied
})

test_that("colour follows truth when truth is supplied, verdict otherwise", {
  pz <- make_pieces()

  # No faulty: colours are the chart's own verdict.
  p_dec <- plot_method_comparison(pz$mon_r, pz$mon_c, pz$ucl_r, pz$ucl_c)
  expect_equal(levels(points_layer(p_dec)$Status), c("No signal", "Signalled"))

  # With faulty: colours are ground truth, as in Figure 6 of the paper.
  p_tru <- plot_method_comparison(pz$mon_r, pz$mon_c, pz$ucl_r, pz$ucl_c,
                                  faulty = c("B3", "B5"))
  d <- points_layer(p_tru)
  expect_equal(levels(d$Status), c("In-control batch", "Faulty batch"))
  expect_setequal(unique(d$Batch[d$Status == "Faulty batch"]), c("B3", "B5"))
  # B5 is faulty and the robust panel does not signal it: red below the line.
  robust_B5 <- d[d$Method == "AFM-MCD (robust)" & d$Batch == "B5", ]
  expect_equal(as.character(robust_B5$Status), "Faulty batch")
  expect_false(robust_B5$flagged)

  # colour_by overrides the automatic choice in both directions.
  p_forced <- plot_method_comparison(pz$mon_r, pz$mon_c, pz$ucl_r, pz$ucl_c,
                                     faulty = c("B3"), colour_by = "decision")
  expect_equal(levels(points_layer(p_forced)$Status),
               c("No signal", "Signalled"))

  # Asking for truth without truth is refused, with the reason.
  expect_error(plot_method_comparison(pz$mon_r, pz$mon_c, pz$ucl_r, pz$ucl_c,
                                      colour_by = "truth"),
               "needs the truth")
  expect_error(plot_method_comparison(pz$mon_r, pz$mon_c, pz$ucl_r, pz$ucl_c,
                                      colour_by = "truth"),
               "never in production")
})

test_that("counts speak of detections only when the truth was supplied", {
  pz <- make_pieces()
  p  <- plot_method_comparison(pz$mon_r, pz$mon_c, pz$ucl_r, pz$ucl_c,
                               faulty = c("B3", "B5"), diagnostics = TRUE)
  counts <- unlist(lapply(p$layers, function(l) l$data$label))

  expect_true(any(grepl("detected 1 of 2 faulty", counts)))
  expect_true(any(grepl("false alarm", counts)))

  # A typo in 'faulty' is caught rather than silently ignored.
  expect_error(plot_method_comparison(pz$mon_r, pz$mon_c, pz$ucl_r, pz$ucl_c,
                                      faulty = "B99"),
               "not among the monitored batches")
})

test_that("the caption exists only when the scale needs defending", {
  pz <- make_pieces()

  p_ratio <- plot_method_comparison(pz$mon_r, pz$mon_c, pz$ucl_r, pz$ucl_c)
  expect_null(p_ratio$labels$caption)
  expect_equal(points_layer(p_ratio)$y,
               points_layer(p_ratio)$T2 / points_layer(p_ratio)$UCL)

  p_t2 <- plot_method_comparison(pz$mon_r, pz$mon_c, pz$ucl_r, pz$ucl_c,
                                 scale = "T2")
  expect_equal(points_layer(p_t2)$y, points_layer(p_t2)$T2)
  expect_true(grepl("NOT comparable", p_t2$labels$caption))
})

test_that("the near-miss band is off by default and labelled when on", {
  pz <- make_pieces()

  p0 <- plot_method_comparison(pz$mon_r, pz$mon_c, pz$ucl_r, pz$ucl_c)
  expect_equal(n_geom(p0, "GeomRect"), 0L)

  p <- plot_method_comparison(pz$mon_r, pz$mon_c, pz$ucl_r, pz$ucl_c,
                              near_miss = 0.20)
  expect_equal(n_geom(p, "GeomRect"), 1L)

  band <- p$layers[[1]]$data
  expect_equal(as.character(band$Method), "Classical Hotelling")
  expect_equal(band$ymin, 0.8)
  expect_equal(band$ymax, 1)
  # Its definition travels with the band, inline, not in a caption.
  expect_equal(band$label, "within 20% of the limit")
  expect_null(p$labels$caption)
})

test_that("x-axis labels stay real batch identifiers, thinned when many", {
  pz <- make_pieces()
  x5 <- plot_method_comparison(pz$mon_r, pz$mon_c,
                               pz$ucl_r, pz$ucl_c)$scales$get_scales("x")
  expect_equal(x5$labels, paste0("B", 1:5))   # few batches: all of them

  study <- make_study()                        # 8 batches, still all
  x8 <- plot_method_comparison(study)$scales$get_scales("x")
  expect_equal(length(x8$labels), 8L)

  # Twenty batches: thinned to every second one, and still identifiers.
  mon_r <- data.frame(Batch = sprintf("F2_B%02d", 1:20), T2 = seq_len(20))
  mon_c <- data.frame(Batch = sprintf("F2_B%02d", 1:20), T2 = seq_len(20))
  x20 <- plot_method_comparison(mon_r, mon_c, 10, 10)$scales$get_scales("x")
  expect_equal(x20$labels, sprintf("F2_B%02d", seq(1, 20, by = 2)))
  expect_false(any(grepl("^[0-9]+$", x20$labels)))
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
                                      diagnostics = "yes"),
               "'diagnostics' must be")
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
