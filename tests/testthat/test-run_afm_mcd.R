make_sim <- function() {
  simulate_batch_process(
    K1 = 12, K2 = 8, I = 20, J = 4,
    outlier_batches_F1 = 1, prop_contam_F1 = 0.07,
    prop_ooc_F2 = 0.30, shift_ooc = 2,
    seed = 20260417
  )
}

test_that("run_afm_mcd returns an afm_mcd_study with every piece in place", {
  sim   <- make_sim()
  vars  <- paste0("Var", 1:4)
  study <- run_afm_mcd(subset(sim, Phase == "Phase 1"),
                       subset(sim, Phase == "Phase 2"),
                       variables = vars, plot = FALSE)

  expect_s3_class(study, "afm_mcd_study")
  expect_named(study, c("calibration", "ucl", "monitoring", "chart",
                        "weight_plot", "classical", "variables", "call"))
  expect_s3_class(study$chart, "ggplot")
  expect_s3_class(study$weight_plot, "ggplot")
  expect_null(study$classical)
  expect_equal(study$variables, vars)
  expect_true("is_ooc" %in% names(study$monitoring))
})

test_that("run_afm_mcd recomputes nothing: the pieces match the direct calls", {
  sim  <- make_sim()
  vars <- paste0("Var", 1:4)
  p1   <- subset(sim, Phase == "Phase 1")
  p2   <- subset(sim, Phase == "Phase 2")

  cal <- calibrate_afm_mcd(p1, vars)
  u   <- ucl_F_adjusted(cal, I = cal$I_phase1)
  mon <- monitor_afm_mcd(p2, cal, vars, ucl = u$UCL)

  study <- run_afm_mcd(p1, p2, variables = vars, plot = FALSE)

  expect_identical(study$calibration, cal)
  expect_identical(study$ucl, u)
  expect_identical(study$monitoring, mon)
})

test_that("run_afm_mcd takes I from the calibration, not from the caller", {
  sim  <- make_sim()
  vars <- paste0("Var", 1:4)
  study <- run_afm_mcd(subset(sim, Phase == "Phase 1"),
                       subset(sim, Phase == "Phase 2"),
                       variables = vars, plot = FALSE)

  expect_equal(study$ucl$parameters$I, study$calibration$I_phase1)
  expect_equal(study$ucl$parameters$I, 20L)
})

test_that("auto-detection always announces the variables it chose", {
  sim  <- make_sim()
  p1   <- subset(sim, Phase == "Phase 1")
  p2   <- subset(sim, Phase == "Phase 2")

  # The message fires even though nobody prints the object.
  expect_message(run_afm_mcd(p1, p2, plot = FALSE),
                 "Auto-detected 4 process variable\\(s\\): Var1, Var2, Var3, Var4")
  # And it warns about exactly the failure mode that motivates it.
  expect_message(run_afm_mcd(p1, p2, plot = FALSE), "identifier or index")

  # Naming the variables explicitly is silent.
  expect_silent(run_afm_mcd(p1, p2, variables = paste0("Var", 1:4),
                            plot = FALSE))
})

test_that("auto-detection would pick up numeric identifier columns", {
  # Tennessee Eastman shape: numeric columns that are not measurements.
  # A *varying* index is the dangerous case, because nothing downstream
  # objects to it: the study runs to completion and looks credible.
  sim <- make_sim()
  p1  <- subset(sim, Phase == "Phase 1")
  p2  <- subset(sim, Phase == "Phase 2")
  p1$sampleIndex <- seq_len(nrow(p1))
  p2$sampleIndex <- seq_len(nrow(p2))

  expect_message(study <- run_afm_mcd(p1, p2, plot = FALSE), "sampleIndex")
  expect_true("sampleIndex" %in% study$variables)
  expect_equal(length(study$variables), 5L)
  # It really does produce a finished, plausible-looking result. The
  # message is the only thing standing between the user and this.
  expect_true(all(is.finite(study$monitoring$T2)))

  # A constant identifier is caught instead: Sw becomes singular. covMcd
  # warns once per batch on the way, hence the suppressWarnings().
  p1$runId <- 1L
  p2$runId <- 1L
  expect_error(
    suppressWarnings(suppressMessages(run_afm_mcd(p1, p2, plot = FALSE))),
    "singular and cannot be inverted"
  )
})

test_that("run_afm_mcd adds the classical baseline only when asked", {
  sim  <- make_sim()
  vars <- paste0("Var", 1:4)
  p1   <- subset(sim, Phase == "Phase 1")
  p2   <- subset(sim, Phase == "Phase 2")

  plain <- run_afm_mcd(p1, p2, variables = vars, plot = FALSE)
  expect_null(plain$classical)

  both <- run_afm_mcd(p1, p2, variables = vars, plot = FALSE,
                      compare_classical = TRUE)
  expect_named(both$classical, c("calibration", "ucl", "monitoring"))
  expect_true("is_ooc" %in% names(both$classical$monitoring))

  # The robust results are untouched by asking for the comparison.
  expect_identical(both$monitoring, plain$monitoring)
  expect_identical(both$ucl, plain$ucl)
})

test_that("print and summary describe the study without altering it", {
  sim   <- make_sim()
  vars  <- paste0("Var", 1:4)
  study <- run_afm_mcd(subset(sim, Phase == "Phase 1"),
                       subset(sim, Phase == "Phase 2"),
                       variables = vars, plot = FALSE)

  expect_output(print(study), "AFM-MCD robust T2 study")
  expect_output(print(study), "Control limit:   UCL =")
  expect_output(print(study), "out of control")

  s <- summary(study)
  expect_s3_class(s, "summary.afm_mcd_study")
  expect_equal(s$n_ooc, sum(study$monitoring$is_ooc))
  expect_equal(s$K1, length(study$calibration$weights))
  expect_equal(s$uniform, 1 / s$K1)

  # The OOC table is sorted by decreasing T2 and carries the ratio.
  if (s$n_ooc > 1) {
    expect_false(is.unsorted(rev(s$ooc_table$T2)))
    expect_equal(s$ooc_table$ratio, s$ooc_table$T2 / study$ucl$UCL)
  }

  expect_output(print(s), "PHASE 1 CALIBRATION")
  expect_output(print(s), "CONTROL LIMIT")
  expect_output(print(s), "NEXT STEP")

  # The summary must not promise a tolerable percentage of bad points, and
  # must not present the weights as a contamination verdict.
  out <- paste(capture.output(print(s)), collapse = " ")
  expect_false(grepl("can be bad", out, fixed = TRUE))
  expect_true(grepl("not contamination", out, fixed = TRUE))

  # No classical block unless it was requested.
  expect_false(grepl("CLASSICAL BASELINE", out, fixed = TRUE))
})

test_that("summary prints the classical block with no verdict attached", {
  sim  <- make_sim()
  vars <- paste0("Var", 1:4)
  s <- summary(run_afm_mcd(subset(sim, Phase == "Phase 1"),
                           subset(sim, Phase == "Phase 2"),
                           variables = vars, plot = FALSE,
                           compare_classical = TRUE))
  out <- paste(capture.output(print(s)), collapse = " ")

  expect_true(grepl("CLASSICAL BASELINE", out, fixed = TRUE))
  expect_true(grepl("trace(Sp) / trace(Sw)", out, fixed = TRUE))
  expect_true(grepl("det ratio", out, fixed = TRUE))
  expect_true(grepl("reference centres", out, fixed = TRUE))
  # Numbers only: no claim that either method is the better one.
  expect_false(grepl("better", out, ignore.case = TRUE))
  expect_false(grepl("outperform", out, ignore.case = TRUE))
})

test_that("the classical block reports three independent measures", {
  sim  <- make_sim()
  vars <- paste0("Var", 1:4)
  study <- run_afm_mcd(subset(sim, Phase == "Phase 1"),
                       subset(sim, Phase == "Phase 2"),
                       variables = vars, plot = FALSE,
                       compare_classical = TRUE)
  s <- summary(study)

  Sp <- study$classical$calibration$Sp
  Sw <- study$calibration$Sw

  # Each measure is what it claims to be, computed from the two matrices.
  expect_equal(s$classical$det_ratio, det(Sp) / det(Sw))
  expect_equal(s$classical$trace_ratio, sum(diag(Sp)) / sum(diag(Sw)))
  expect_equal(s$classical$centre_gap,
               sqrt(sum((study$classical$calibration$mu_global -
                           study$calibration$mu_r)^2)))

  # The three are independent: the determinant ratio can stay near 1 while
  # the trace ratio does not, because contamination inflates variances and
  # correlations at once and the two effects cancel in the determinant.
  expect_true(all(is.finite(c(s$classical$det_ratio,
                              s$classical$trace_ratio,
                              s$classical$centre_gap))))
})

test_that("run_afm_mcd writes figures only when save_path is given", {
  sim  <- make_sim()
  vars <- paste0("Var", 1:4)
  stem <- file.path(tempdir(), "test_study")
  files <- paste0(stem, c("_control_chart.png", "_control_chart.pdf",
                          "_afm_weights.png", "_afm_weights.pdf"))
  on.exit(unlink(files), add = TRUE)

  run_afm_mcd(subset(sim, Phase == "Phase 1"),
              subset(sim, Phase == "Phase 2"),
              variables = vars, plot = FALSE)
  expect_false(any(file.exists(files)))

  run_afm_mcd(subset(sim, Phase == "Phase 1"),
              subset(sim, Phase == "Phase 2"),
              variables = vars, plot = FALSE, save_path = stem)
  expect_true(all(file.exists(files)))
})

test_that("run_afm_mcd validates its inputs", {
  sim  <- make_sim()
  vars <- paste0("Var", 1:4)
  p1   <- subset(sim, Phase == "Phase 1")
  p2   <- subset(sim, Phase == "Phase 2")

  expect_error(run_afm_mcd("nope", p2, vars), "'phase1' must be a data frame")
  expect_error(run_afm_mcd(p1, "nope", vars), "'phase2' must be a data frame")
  expect_error(run_afm_mcd(p1[setdiff(names(p1), "Batch")], p2, vars),
               "Column 'Batch' not found in 'phase1'")
  expect_error(run_afm_mcd(p1, p2, "Missing"), "not found in 'phase1'")
  expect_error(run_afm_mcd(p1, p2[setdiff(names(p2), "Var4")],
                           vars),
               "not found in 'phase2'")
  expect_error(run_afm_mcd(p1, p2, vars, plot = "yes"), "'plot' must be")
  expect_error(run_afm_mcd(p1, p2, vars, compare_classical = 1),
               "'compare_classical' must be")
  expect_error(run_afm_mcd(p1, p2, vars, save_path = ""), "save_path")
})
