# Naming the batch column "Batch" is this package's default, not a
# requirement of the method. These tests cover the five functions that read a
# user-supplied data frame.

renamed <- function(d, new = "Lote") {
  names(d)[names(d) == "Batch"] <- new
  d
}

test_that("batch_col reads a differently named column, in all five functions", {
  p1 <- renamed(afm_phase1)
  p2 <- renamed(afm_phase2)
  vars <- paste0("Var", 1:4)

  cal_ref <- calibrate_afm_mcd(afm_phase1, vars)
  cal_new <- calibrate_afm_mcd(p1, vars, batch_col = "Lote")
  expect_equal(cal_new$mu_r, cal_ref$mu_r)
  expect_equal(cal_new$Sw, cal_ref$Sw)
  expect_equal(cal_new$weights, cal_ref$weights)

  mon_ref <- monitor_afm_mcd(afm_phase2, cal_ref, vars)
  mon_new <- monitor_afm_mcd(p2, cal_new, vars, batch_col = "Lote")
  expect_equal(mon_new$T2, mon_ref$T2)

  cc_ref <- hotelling_classical_calibrate(afm_phase1, vars)
  cc_new <- hotelling_classical_calibrate(p1, vars, batch_col = "Lote")
  expect_equal(cc_new$Sp, cc_ref$Sp)
  expect_equal(cc_new$mu_global, cc_ref$mu_global)

  cm_ref <- hotelling_classical_monitor(afm_phase2, cc_ref, vars)
  cm_new <- hotelling_classical_monitor(p2, cc_new, vars, batch_col = "Lote")
  expect_equal(cm_new$T2, cm_ref$T2)

  st_ref <- run_afm_mcd(afm_phase1, afm_phase2, vars, plot = FALSE)
  st_new <- run_afm_mcd(p1, p2, vars, plot = FALSE, batch_col = "Lote")
  expect_equal(st_new$monitoring$T2, st_ref$monitoring$T2)
  expect_equal(st_new$ucl$UCL, st_ref$ucl$UCL)
})

test_that("the output column stays 'Batch' whatever the input was called", {
  # Renaming what is read must not rename what is returned: the plots and the
  # study object depend on that name, and it is frozen output.
  p1 <- renamed(afm_phase1)
  p2 <- renamed(afm_phase2)
  vars <- paste0("Var", 1:4)

  cal <- calibrate_afm_mcd(p1, vars, batch_col = "Lote")
  mon <- monitor_afm_mcd(p2, cal, vars, batch_col = "Lote")
  expect_named(mon, c("Batch", "I", "T2"))

  cc  <- hotelling_classical_calibrate(p1, vars, batch_col = "Lote")
  cm  <- hotelling_classical_monitor(p2, cc, vars, batch_col = "Lote")
  expect_named(cm, c("Batch", "I", "T2"))

  # And the identifiers themselves survive intact.
  expect_equal(mon$Batch, as.character(unique(afm_phase2$Batch)))

  # So the plotting functions keep working with a renamed input.
  st <- run_afm_mcd(p1, p2, vars, plot = FALSE, compare_classical = TRUE,
                    batch_col = "Lote")
  expect_s3_class(plot_control_chart(st$monitoring, UCL = st$ucl$UCL), "ggplot")
  expect_s3_class(plot_method_comparison(st), "ggplot")
})

test_that("a missing batch column says what it looked for and what is there", {
  p1 <- renamed(afm_phase1)
  vars <- paste0("Var", 1:4)

  err <- tryCatch(calibrate_afm_mcd(p1, vars), error = function(e) e)
  msg <- conditionMessage(err)
  expect_true(grepl("Column 'Batch' not found in 'data'", msg, fixed = TRUE))
  expect_true(grepl("The available columns are: Lote, Phase, Status", msg,
                    fixed = TRUE))
  expect_true(grepl("batch_col = \"<name>\"", msg, fixed = TRUE))

  # Every entry point names itself in the hint, not some other function.
  expect_error(monitor_afm_mcd(renamed(afm_phase2),
                               calibrate_afm_mcd(afm_phase1, vars), vars),
               "monitor_afm_mcd\\(new_data, calibration, variables, batch_col")
  expect_error(hotelling_classical_calibrate(p1, vars),
               "hotelling_classical_calibrate\\(data, variables, batch_col")
  expect_error(run_afm_mcd(p1, renamed(afm_phase2), vars),
               "run_afm_mcd\\(phase1, phase2, batch_col")

  # And a bad batch_col is rejected before anything else happens.
  expect_error(calibrate_afm_mcd(afm_phase1, vars, batch_col = c("a", "b")),
               "single non-empty column name")
  expect_error(calibrate_afm_mcd(afm_phase1, vars, batch_col = ""),
               "single non-empty column name")
})

test_that("auto-detection excludes the renamed batch column", {
  p1 <- renamed(afm_phase1)
  p2 <- renamed(afm_phase2)

  expect_message(
    study <- run_afm_mcd(p1, p2, plot = FALSE, batch_col = "Lote"),
    "Auto-detected 4 process variable"
  )
  expect_equal(study$variables, paste0("Var", 1:4))
  expect_false("Lote" %in% study$variables)
})
