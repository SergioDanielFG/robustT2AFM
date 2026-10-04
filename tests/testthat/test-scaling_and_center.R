# Tests of the two options added in version 0.3.0 (review round 1):
# scaling = "mad" (comment 4) and center = "mcd" (comments 6 and 7).

vars <- paste0("Var", 1:4)

test_that("scaling = 'none', center = 'mean' reproduces version 0.2.0", {
  set.seed(1)
  cal <- calibrate_afm_mcd(afm_phase1, vars, scaling = "none", center = "mean")

  # 0.2.0: raw first eigenvalues, inverse weights, simple average of centers.
  lam <- sapply(cal$mcd_covariances, function(S)
    max(eigen(S, symmetric = TRUE, only.values = TRUE)$values))
  expect_equal(unname(cal$lambda1), unname(lam))
  expect_equal(unname(cal$weights), unname((1 / lam) / sum(1 / lam)))
  expect_equal(unname(cal$mu_r),
               unname(colMeans(do.call(rbind, cal$mcd_centers))))
  expect_equal(unname(cal$scale), rep(1, 4))
})

test_that("with scaling = 'mad' a change of units leaves weights and T2 unchanged", {
  set.seed(2)
  cal1 <- calibrate_afm_mcd(afm_phase1, vars)
  ph1 <- afm_phase1; ph1$Var1 <- 100 * ph1$Var1
  ph2 <- afm_phase2; ph2$Var1 <- 100 * ph2$Var1
  set.seed(2)
  cal2 <- calibrate_afm_mcd(ph1, vars)

  expect_equal(cal2$weights, cal1$weights, tolerance = 1e-10)
  expect_equal(cal2$mu_r[["Var1"]], 100 * cal1$mu_r[["Var1"]], tolerance = 1e-10)
  expect_equal(monitor_afm_mcd(ph2, cal2, vars)$T2,
               monitor_afm_mcd(afm_phase2, cal1, vars)$T2, tolerance = 1e-8)
})

test_that("without scaling the weights do depend on the units", {
  set.seed(3)
  cal1 <- calibrate_afm_mcd(afm_phase1, vars, scaling = "none")
  ph1 <- afm_phase1; ph1$Var1 <- 100 * ph1$Var1
  set.seed(3)
  cal2 <- calibrate_afm_mcd(ph1, vars, scaling = "none")
  expect_false(isTRUE(all.equal(cal2$weights, cal1$weights)))
})

test_that("center = 'mcd' resists whole batches shifted in mean", {
  ph1 <- afm_phase1
  moved <- unique(ph1$Batch)[1:6]
  ph1[ph1$Batch %in% moved, vars] <- ph1[ph1$Batch %in% moved, vars] + 5

  set.seed(4)
  cal_mcd  <- calibrate_afm_mcd(ph1, vars, center = "mcd")
  set.seed(4)
  cal_mean <- calibrate_afm_mcd(ph1, vars, center = "mean")

  # 6 of 30 batches moved 5 units: the average moves about 1, the MCD barely.
  expect_gt(max(abs(cal_mean$mu_r)), 0.7)
  expect_lt(max(abs(cal_mcd$mu_r)), 0.3)
})

test_that("invalid options stop with a clear message", {
  expect_error(calibrate_afm_mcd(afm_phase1, vars, scaling = "zscore"))
  expect_error(calibrate_afm_mcd(afm_phase1, vars, center_alpha = 0.3),
               "center_alpha")

  flat <- afm_phase1; flat$Var2 <- 1
  expect_error(calibrate_afm_mcd(flat, vars), "median absolute deviation")

  few <- afm_phase1[afm_phase1$Batch %in% unique(afm_phase1$Batch)[1:4], ]
  expect_error(calibrate_afm_mcd(few, vars, center = "mcd"),
               "at least twice as many valid batches")
  expect_silent(calibrate_afm_mcd(few, vars, center = "mean"))
})

test_that("the MCD of the centers is reproducible and does not disturb the user's RNG", {
  # Same Phase 1 twice, no set.seed(): identical center (internal seed).
  c1 <- calibrate_afm_mcd(afm_phase1, vars)$mu_r
  c2 <- calibrate_afm_mcd(afm_phase1, vars)$mu_r
  expect_identical(c1, c2)

  # The internal seed is undone: after the calibration, the user's random
  # stream is exactly where center = "mean" (no internal seed) leaves it.
  set.seed(99); invisible(calibrate_afm_mcd(afm_phase1, vars, center = "mean"))
  a <- runif(3)
  set.seed(99); invisible(calibrate_afm_mcd(afm_phase1, vars, center = "mcd"))
  b <- runif(3)
  expect_identical(a, b)
})
