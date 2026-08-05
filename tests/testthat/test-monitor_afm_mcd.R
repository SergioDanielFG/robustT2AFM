test_that("monitor_afm_mcd returns non-negative T^2 for every batch", {
  set.seed(20260417)
  sim  <- simulate_batch_process(K1 = 12, K2 = 5, I = 20, J = 4,
                                 prop_ooc_F2 = 0.4, shift_ooc = 2,
                                 seed = 20260417)
  vars <- paste0("Var", 1:4)
  cal  <- suppressMessages(
    calibrate_afm_mcd(subset(sim, Phase == "Phase 1"), vars)
  )
  mon  <- monitor_afm_mcd(subset(sim, Phase == "Phase 2"), cal, vars)

  expect_s3_class(mon, "data.frame")
  expect_named(mon, c("Batch", "I", "T2"))
  expect_true(all(mon$T2 >= 0))
  expect_true(all(mon$I == 20L))
})

test_that("monitor_afm_mcd yields T^2 = 0 when batch mean equals mu_r", {
  # Trivial calibration: mu_r = 0, Sw = I; batch with zero mean must return 0.
  cal <- list(
    mu_r    = c(0, 0, 0, 0),
    Sw      = diag(4),
    weights = setNames(rep(1 / 30, 30), paste0("B", 1:30))
  )
  new_batch <- data.frame(
    Batch = "F2_new",
    Var1  = c(-1, 1),
    Var2  = c(-1, 1),
    Var3  = c(-1, 1),
    Var4  = c(-1, 1)
  )
  mon <- monitor_afm_mcd(new_batch, cal, paste0("Var", 1:4))
  expect_equal(mon$T2, 0)
  expect_equal(mon$I, 2L)
})

test_that("monitor_afm_mcd adds is_ooc only when ucl is supplied", {
  sim  <- simulate_batch_process(K1 = 12, K2 = 5, I = 20, J = 4,
                                 prop_ooc_F2 = 0.4, shift_ooc = 2,
                                 seed = 20260417)
  vars <- paste0("Var", 1:4)
  cal  <- calibrate_afm_mcd(subset(sim, Phase == "Phase 1"), vars)
  ph2  <- subset(sim, Phase == "Phase 2")

  mon_plain <- monitor_afm_mcd(ph2, cal, vars)
  mon_flag  <- monitor_afm_mcd(ph2, cal, vars,
                               ucl = ucl_F_adjusted(cal, I = 20)$UCL)

  # Default output is untouched: three columns, not one more.
  expect_named(mon_plain, c("Batch", "I", "T2"))

  # ucl = NULL reproduces the default output exactly.
  expect_identical(monitor_afm_mcd(ph2, cal, vars, ucl = NULL), mon_plain)

  # The flag is additive: the first three columns are unchanged.
  expect_named(mon_flag, c("Batch", "I", "T2", "is_ooc"))
  expect_identical(mon_flag[c("Batch", "I", "T2")], mon_plain)
  expect_type(mon_flag$is_ooc, "logical")
})

test_that("monitor_afm_mcd flags strictly above the limit, never on it", {
  # mu_r = 0, Sw = I, batch of 2 observations with x_bar = (2, 2, 2, 2):
  # T2 = I * ||x_bar||^2 = 2 * 16 = 32, an exact value with no rounding.
  cal  <- list(mu_r = c(0, 0, 0, 0), Sw = diag(4))
  vars <- paste0("Var", 1:4)
  new_batch <- data.frame(
    Batch = "F2_new",
    Var1  = c(1, 3), Var2 = c(1, 3),
    Var3  = c(1, 3), Var4 = c(1, 3)
  )
  expect_equal(monitor_afm_mcd(new_batch, cal, vars)$T2, 32)

  expect_false(monitor_afm_mcd(new_batch, cal, vars, ucl = 33)$is_ooc)
  expect_true(monitor_afm_mcd(new_batch, cal, vars, ucl = 31)$is_ooc)

  # Exactly on the limit: NOT flagged (strictly greater, not >=).
  expect_false(monitor_afm_mcd(new_batch, cal, vars, ucl = 32)$is_ooc)
})

test_that("monitor_afm_mcd validates ucl", {
  cal  <- list(mu_r = c(0, 0, 0, 0), Sw = diag(4))
  vars <- paste0("Var", 1:4)
  new_batch <- data.frame(
    Batch = "F2_new",
    Var1  = c(-1, 1), Var2 = c(-1, 1),
    Var3  = c(-1, 1), Var4 = c(-1, 1)
  )

  expect_error(monitor_afm_mcd(new_batch, cal, vars, ucl = -1),
               "single positive, finite number")
  expect_error(monitor_afm_mcd(new_batch, cal, vars, ucl = 0),
               "single positive, finite number")
  expect_error(monitor_afm_mcd(new_batch, cal, vars, ucl = c(10, 20)),
               "single positive, finite number")
  expect_error(monitor_afm_mcd(new_batch, cal, vars, ucl = "19.7"),
               "single positive, finite number")
  expect_error(monitor_afm_mcd(new_batch, cal, vars, ucl = NA_real_),
               "single positive, finite number")

  # Passing the whole ucl_F_adjusted() object is the predictable mistake:
  # the error must name the fix.
  expect_error(monitor_afm_mcd(new_batch, cal, vars,
                               ucl = list(UCL = 19.7, method = "F-adjusted")),
               "ucl_F_adjusted\\(calibration, I\\)\\$UCL")
})

test_that("monitor_afm_mcd reports a singular Sw with a usable message", {
  # Var4 is an exact copy of Var3, so Sw is rank-deficient and not invertible.
  Sw <- matrix(c(1, 0, 0, 0,
                 0, 1, 0, 0,
                 0, 0, 1, 1,
                 0, 0, 1, 1), nrow = 4, byrow = TRUE)
  cal <- list(mu_r = c(0, 0, 0, 0), Sw = Sw)
  new_batch <- data.frame(
    Batch = "F2_new",
    Var1  = c(-1, 1), Var2 = c(-1, 1),
    Var3  = c(-1, 1), Var4 = c(-1, 1)
  )

  expect_error(monitor_afm_mcd(new_batch, cal, paste0("Var", 1:4)),
               "singular and cannot be inverted")
  # The message must tell the engineer what to do, not only what failed.
  expect_error(monitor_afm_mcd(new_batch, cal, paste0("Var", 1:4)),
               "redundant")
})
