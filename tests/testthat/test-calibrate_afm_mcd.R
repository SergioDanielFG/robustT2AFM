test_that("calibrate_afm_mcd returns a coherent AFM/MCD reference", {
  set.seed(20260417)
  sim  <- simulate_batch_process(K1 = 12, K2 = 0, I = 20, J = 4,
                                 seed = 20260417)
  vars <- paste0("Var", 1:4)
  cal  <- suppressMessages(calibrate_afm_mcd(sim, vars))

  # Return structure ---------------------------------------------------------
  expect_type(cal, "list")
  expect_named(cal, c("mu_r", "Sw", "weights", "mcd_centers",
                      "mcd_covariances", "lambda1", "mcd_alpha", "I_phase1"))
  # mcd_clean_obs was removed together with the bootstrap UCL.
  expect_false("mcd_clean_obs" %in% names(cal))

  # I_phase1 records the Phase 1 batch size (first valid batch).
  expect_equal(cal$I_phase1, 20L)

  # Default mcd_alpha matches Ruiz-Barzola et al. (2026).
  expect_equal(cal$mcd_alpha, 0.67)

  # AFM inverse weights: strictly positive and sum to exactly 1.
  expect_true(all(cal$weights > 0))
  expect_equal(sum(cal$weights), 1)

  # Sw is J x J, symmetric and positive-definite.
  expect_equal(dim(cal$Sw), c(4L, 4L))
  expect_equal(cal$Sw, t(cal$Sw))
  min_eig <- min(eigen(cal$Sw, symmetric = TRUE,
                       only.values = TRUE)$values)
  expect_gt(min_eig, 0)

  # mu_r has length J.
  expect_length(cal$mu_r, 4L)
})

test_that("calibrate_afm_mcd validates its inputs", {
  df <- data.frame(Batch = "B1", Var1 = 1, Var2 = 2)
  expect_error(calibrate_afm_mcd("not a df", c("Var1", "Var2")),
               "must be a data frame")
  expect_error(calibrate_afm_mcd(data.frame(Var1 = 1), "Var1"),
               "must contain a column named 'Batch'")
  expect_error(calibrate_afm_mcd(df, "Missing"), "not found in data")
  expect_error(calibrate_afm_mcd(df, c("Var1", "Var2"), mcd_alpha = 0.5),
               "mcd_alpha must be")
})
