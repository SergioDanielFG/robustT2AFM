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
