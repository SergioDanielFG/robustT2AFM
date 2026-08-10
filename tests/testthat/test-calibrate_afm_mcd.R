test_that("calibrate_afm_mcd returns a coherent AFM/MCD reference", {
  set.seed(20260417)
  sim  <- simulate_batch_process(K1 = 12, K2 = 0, I = 20, J = 4,
                                 seed = 20260417)
  vars <- paste0("Var", 1:4)
  cal  <- suppressMessages(calibrate_afm_mcd(sim, vars))

  # Return structure ---------------------------------------------------------
  expect_type(cal, "list")
  expect_named(cal, c("mu_r", "Sw", "weights", "mcd_centers",
                      "mcd_covariances", "lambda1", "mcd_alpha", "I_phase1",
                      "batch_sizes"))
  # mcd_clean_obs was removed together with the bootstrap UCL.
  expect_false("mcd_clean_obs" %in% names(cal))

  # I_phase1 records the Phase 1 batch size (first valid batch).
  expect_equal(cal$I_phase1, 20L)

  # Default mcd_alpha matches Frutos-Galarza et al. (2026).
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

test_that("unequal Phase 1 batches are announced and not silently averaged over", {
  sim  <- simulate_batch_process(K1 = 6, K2 = 0, I = 20, J = 4,
                                 seed = 20260425)
  vars <- paste0("Var", 1:4)

  # Equal batches: silent, and I_phase1 is that size.
  cal_eq <- calibrate_afm_mcd(sim, vars)
  expect_equal(cal_eq$I_phase1, 20L)
  expect_equal(unname(cal_eq$batch_sizes), rep(20L, 6))

  # Two shortened batches, as happens when a run stops early or rows are
  # discarded. Sizes 14, 17, 20, 20, 20, 20: the per-batch m* are 9, 11 and
  # 13 x 4, whose mean is 12, and the I giving m* = 12 is round(12/0.67) = 18.
  drop <- c(rownames(sim[sim$Batch == "F1_B01", ])[1:6],
            rownames(sim[sim$Batch == "F1_B02", ])[1:3])
  uneq <- sim[!rownames(sim) %in% drop, ]

  expect_warning(cal_un <- calibrate_afm_mcd(uneq, vars),
                 "do not all have the same number of observations")
  expect_warning(calibrate_afm_mcd(uneq, vars), "range from 14 to 20")
  expect_warning(calibrate_afm_mcd(uneq, vars), "4 of the 6 batches have 20")
  expect_warning(calibrate_afm_mcd(uneq, vars), "equals 12")
  expect_warning(calibrate_afm_mcd(uneq, vars), "ucl_F_adjusted\\(calibration")

  expect_equal(unname(cal_un$batch_sizes), c(14L, 17L, 20L, 20L, 20L, 20L))
  expect_equal(cal_un$I_phase1, 18L)          # not the first batch (14)

  # The recorded size no longer depends on which batch happens to come first.
  rev_order <- uneq[order(uneq$Batch, decreasing = TRUE), ]
  expect_equal(suppressWarnings(
    calibrate_afm_mcd(rev_order, vars)$I_phase1), cal_un$I_phase1)
})

test_that("I_phase1 follows mean(m*_k), not the rounded mean batch size", {
  sim  <- simulate_batch_process(K1 = 6, K2 = 0, I = 20, J = 4,
                                 seed = 20260425)
  vars <- paste0("Var", 1:4)

  # Sizes 11, 11, 12, 20, 20, 20 is a case where the two rules disagree:
  #   mean(round(I_k * h))       = (7 + 7 + 8 + 13*3)/6 = 10.167 -> m* = 10
  #     -> I = round(10/0.67)    = 15
  #   round(mean(I_k) * h)       = round(round(15.667) * 0.67) = 11  -> I = 16
  drop <- c(rownames(sim[sim$Batch == "F1_B01", ])[1:9],
            rownames(sim[sim$Batch == "F1_B02", ])[1:9],
            rownames(sim[sim$Batch == "F1_B03", ])[1:8])
  uneq <- sim[!rownames(sim) %in% drop, ]

  cal <- suppressWarnings(calibrate_afm_mcd(uneq, vars))
  expect_equal(unname(cal$batch_sizes), c(11L, 11L, 12L, 20L, 20L, 20L))
  expect_equal(cal$I_phase1, 15L)             # the m*-matching size
  expect_false(cal$I_phase1 == 16L)           # not the rounded mean size

  # And the recorded size really does reproduce the intended m*.
  expect_equal(round(cal$I_phase1 * cal$mcd_alpha), 10)
  expect_equal(suppressWarnings(
    ucl_F_adjusted(cal, I = cal$I_phase1)$parameters$m_star), 10)
})

test_that("ucl_F_adjusted reports that an unequal calibration has no exact I", {
  cal <- list(Sw = diag(4), mcd_alpha = 0.67,
              weights = setNames(rep(1 / 30, 30), paste0("B", 1:30)),
              I_phase1 = 18L,
              batch_sizes = c(rep(20L, 28), 14L, 17L))

  expect_warning(ucl_F_adjusted(cal, I = 18),
                 "not equal \\(14 to 20\\).*approximation")

  # With equal sizes the historical mismatch warning is untouched.
  cal$batch_sizes <- rep(20L, 30)
  cal$I_phase1    <- 20L
  expect_warning(ucl_F_adjusted(cal, I = 15), "does not match the Phase 1")
  expect_silent(ucl_F_adjusted(cal, I = 20))
})

test_that("calibrate_afm_mcd validates its inputs", {
  df <- data.frame(Batch = "B1", Var1 = 1, Var2 = 2)
  expect_error(calibrate_afm_mcd("not a df", c("Var1", "Var2")),
               "must be a data frame")
  expect_error(calibrate_afm_mcd(data.frame(Var1 = 1), "Var1"),
               "Column 'Batch' not found in 'data'")
  # The "not found" error lists what is there, like its batch_col sibling.
  expect_error(calibrate_afm_mcd(df, "Missing"),
               "Variable\\(s\\) not found in 'data': Missing")
  expect_error(calibrate_afm_mcd(df, "Missing"),
               "numeric columns available are: Var1, Var2")
  expect_error(calibrate_afm_mcd(df, c("Var1", "Var2"), mcd_alpha = 0.5),
               "mcd_alpha must be")
  expect_error(calibrate_afm_mcd(df, c("Var1", "Var2"), verbose = "yes"),
               "'verbose' must be a single logical")
})

test_that("calibrate_afm_mcd rejects a single variable and says what to do", {
  vars <- paste0("Var", 1:4)

  expect_error(calibrate_afm_mcd(afm_phase1, "Var1"),
               "needs at least 2 process variables")
  # El mensaje tiene que ofrecer la alternativa, no solo negarse.
  expect_error(calibrate_afm_mcd(afm_phase1, "Var1"),
               "Shewhart X-bar chart")

  # El guardia va despues de los dos que ya existian: una columna de lote
  # ausente o una variable inexistente siguen dando su propio error.
  expect_error(calibrate_afm_mcd(data.frame(Var1 = 1), "Var1"),
               "Column 'Batch' not found in 'data'")
  expect_error(calibrate_afm_mcd(afm_phase1, "Missing"),
               "Variable\\(s\\) not found in 'data'")

  # Dos variables siguen bastando: el limite es J < 2, no J < 4.
  expect_type(calibrate_afm_mcd(afm_phase1, vars[1:2]), "list")
})

test_that("calibrate_afm_mcd rejects a non-numeric variable instead of ranking it", {
  vars <- paste0("Var", 1:4)

  d <- afm_phase1
  d$Var2 <- as.character(d$Var2)
  expect_error(calibrate_afm_mcd(d, vars), "not numeric: Var2")

  # Lo que se esta impidiendo: sin el guardia, data.matrix() convertia la
  # columna de texto en codigos de factor ordenados alfabeticamente y la
  # calibracion devolvia numeros finitos y creibles. El centro de Var2 se
  # iba de ~0.01 a ~10.3 sin que nada avisara.
  cal_ok <- calibrate_afm_mcd(afm_phase1, vars)
  expect_lt(abs(cal_ok$mu_r[["Var2"]]), 1)
})

test_that("calibrate_afm_mcd rejects non-finite values naming the batch", {
  vars <- paste0("Var", 1:4)

  d <- afm_phase1
  d$Var1[5] <- NA
  expect_error(calibrate_afm_mcd(d, vars),
               "Batch 'F1_B01' contains non-finite values")

  d2 <- afm_phase1
  d2$Var3[nrow(d2)] <- Inf
  expect_error(calibrate_afm_mcd(d2, vars), "contains non-finite values")
})

test_that("calibrate_afm_mcd is silent by default and talks under verbose", {
  sim  <- simulate_batch_process(K1 = 5, K2 = 0, I = 20, J = 4,
                                 seed = 20260417)
  vars <- paste0("Var", 1:4)

  # Default: no message at all (scripts no longer need suppressMessages()).
  expect_silent(calibrate_afm_mcd(sim, vars))

  # verbose = TRUE restores the historical listing of valid batches.
  expect_message(calibrate_afm_mcd(sim, vars, verbose = TRUE),
                 "Valid batches used for calibration \\(5\\)")

  # verbose does not touch the numbers.
  cal_quiet <- calibrate_afm_mcd(sim, vars)
  cal_loud  <- suppressMessages(calibrate_afm_mcd(sim, vars, verbose = TRUE))
  expect_identical(cal_quiet, cal_loud)
})
