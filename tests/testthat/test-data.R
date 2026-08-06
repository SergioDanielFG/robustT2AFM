test_that("the shipped data sets are exactly what their @source says", {
  # The call below is the one documented in ?afm_phase1 and in
  # data-raw/make_datasets.R. If this test ever fails on a future R, the
  # likely cause is a change in the random number generator or in
  # MASS::mvrnorm, not a defect in the data: the shipped objects are the
  # ones that reproduce every example and figure in the package, so they
  # are authoritative and this test would be what needs revisiting.
  sim <- simulate_batch_process(
    K1 = 30, K2 = 20, I = 20, J = 4, rho = 0.6,
    outlier_batches_F1 = 6, outlier_rate = 0.20, outlier_shift = 4,
    prop_ooc_F2 = 0.5, shift_ooc = 1.0,
    seed = 20260425
  )
  p1 <- droplevels(subset(sim, Phase == "Phase 1"))
  p2 <- droplevels(subset(sim, Phase == "Phase 2"))
  rownames(p1) <- NULL
  rownames(p2) <- NULL

  expect_identical(p1, afm_phase1)
  expect_identical(p2, afm_phase2)
})

test_that("the shipped data sets have the shape the documentation promises", {
  expect_s3_class(afm_phase1, "data.frame")
  expect_s3_class(afm_phase2, "data.frame")

  cols <- c("Batch", "Phase", "Status", "ContaminationType",
            paste0("Var", 1:4))
  expect_named(afm_phase1, cols)
  expect_named(afm_phase2, cols)

  expect_equal(nrow(afm_phase1), 600L)
  expect_equal(nrow(afm_phase2), 400L)
  expect_equal(nlevels(afm_phase1$Batch), 30L)
  expect_equal(nlevels(afm_phase2$Batch), 20L)
  expect_true(all(table(afm_phase1$Batch) == 20L))
  expect_true(all(table(afm_phase2$Batch) == 20L))

  # 6 contaminated Phase 1 batches, 10 off-target Phase 2 batches.
  contam <- unique(afm_phase1$Batch[afm_phase1$ContaminationType == "Outliers"])
  expect_length(contam, 6L)
  ooc <- unique(afm_phase2$Batch[afm_phase2$Status == "Out of Control"])
  expect_length(ooc, 10L)
})

test_that("the ground-truth columns cannot be mistaken for process variables", {
  # Status, ContaminationType, Phase and Batch are all factors, so the
  # auto-detection in run_afm_mcd picks exactly the four Var columns and the
  # examples need no 'variables' argument.
  for (d in list(afm_phase1, afm_phase2)) {
    numeric_cols <- names(d)[vapply(d, is.numeric, logical(1))]
    expect_equal(numeric_cols, paste0("Var", 1:4))
  }
})
